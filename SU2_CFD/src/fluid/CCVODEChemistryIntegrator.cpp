#include "../../include/fluid/CCVODEChemistryIntegrator.hpp"

#include <algorithm>
#include <limits>
#include <utility>
#include <vector>

#include "../../include/fluid/CCanteraReactiveGas.hpp"
#include "../../../Common/include/CConfig.hpp"

#ifdef HAVE_SUNDIALS
#include "cvode/cvode.h"
#include "cvode/cvode_ls.h"
#include "nvector/nvector_serial.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "sunmatrix/sunmatrix_dense.h"
#include "sundials/sundials_context.h"
#endif

namespace {

#ifdef HAVE_SUNDIALS
struct ChemistryRHSData {
  CCanteraReactiveGas* reactiveGas = nullptr;
  su2double density = 0.0;
  std::vector<su2double> massFractions;
  std::vector<su2double> wdot;
  std::vector<su2double> speciesEnthalpies;
  std::vector<su2double> speciesInternalEnergies;
};

bool CVodeCallSucceeded(int flag) { return flag >= 0; }

int DetonationChemistryRHS(sunrealtype, N_Vector state, N_Vector rhs, void* userData) {
  auto* data = static_cast<ChemistryRHSData*>(userData);
  auto* stateData = N_VGetArrayPointer(state);
  auto* rhsData = N_VGetArrayPointer(rhs);

  const unsigned short nSpecies = static_cast<unsigned short>(data->massFractions.size());
  const su2double temperature = max(static_cast<su2double>(stateData[0]), su2double(1.0));

  for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    data->massFractions[iSpecies] = max(static_cast<su2double>(stateData[iSpecies + 1]), su2double(0.0));
  }
  data->reactiveGas->NormalizeMassFractions(data->massFractions, 0.0);

  su2double cvMix = 0.0;
  su2double qdot = 0.0;
  data->reactiveGas->ComputeChemistrySourceTerms(data->density, temperature, data->massFractions, data->wdot,
                                                 data->speciesEnthalpies, data->speciesInternalEnergies, cvMix, qdot);

  rhsData[0] = 0.0;
  const su2double rhoCv = max(data->density * cvMix, su2double(1.0e-30));
  for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    rhsData[iSpecies + 1] = data->wdot[iSpecies] / data->density;
    rhsData[0] -= data->speciesInternalEnergies[iSpecies] * data->wdot[iSpecies] / rhoCv;
  }

  return 0;
}
#endif

}  // namespace

CCVODEChemistryIntegrator::CCVODEChemistryIntegrator(CCanteraReactiveGas& reactiveGasIn, const CConfig* configIn)
    : reactiveGas(reactiveGasIn), config(configIn) {
#ifndef HAVE_SUNDIALS
  SU2_MPI::Error("DETONATION_EULER was selected but this binary was not compiled with SUNDIALS/CVODE support. Reconfigure with -Denable-sundials=true.", CURRENT_FUNCTION);
#endif
}

ChemistryAdvanceResult CCVODEChemistryIntegrator::IntegrateCell(su2double dt,
                                                                su2double density,
                                                                su2double temperature,
                                                                const std::vector<su2double>& massFractions) {
  ChemistryAdvanceResult result;
  result.density = density;
  result.temperatureBefore = temperature;
  result.temperature = min(max(temperature, config->GetChemistry_Clip_TMin()), config->GetChemistry_Clip_TMax());
  result.temperatureClipCorrection = fabs(result.temperature - temperature);
  result.massFractions = massFractions;
  result.sumYBefore = 0.0;
  result.minYBefore = std::numeric_limits<su2double>::max();
  for (auto value : massFractions) {
    result.sumYBefore += value;
    result.minYBefore = min(result.minYBefore, value);
  }
  if (result.minYBefore == std::numeric_limits<su2double>::max()) result.minYBefore = 0.0;
  reactiveGas.NormalizeMassFractions(result.massFractions, config->GetChemistry_Clip_YMin());

  const auto evaluateState = [&](su2double temperatureState,
                                 const std::vector<su2double>& massFractionsState,
                                 su2double& internalEnergy,
                                 su2double& pressure) {
    reactiveGas.SetMassFractions(massFractionsState.data(), reactiveGas.GetNumSpecies());
    reactiveGas.SetTDState_rhoT(density, temperatureState);
    internalEnergy = reactiveGas.GetStaticEnergy();
    pressure = reactiveGas.GetPressure();
  };

  evaluateState(result.temperature, result.massFractions, result.internalEnergyBefore, result.pressureBefore);
  const su2double eOld = result.internalEnergyBefore;

  if (dt <= 0.0 || density <= 0.0) {
    result.success = true;
    result.substeps = 1;
    result.qdot = 0.0;
    result.sumYAfter = 0.0;
    result.minYAfter = std::numeric_limits<su2double>::max();
    for (auto value : result.massFractions) {
      result.sumYAfter += value;
      result.minYAfter = min(result.minYAfter, value);
    }
    if (result.minYAfter == std::numeric_limits<su2double>::max()) result.minYAfter = 0.0;
    result.pressureAfter = result.pressureBefore;
    result.internalEnergyAfter = result.internalEnergyBefore;
    return result;
  }

#ifndef HAVE_SUNDIALS
  return result;
#else
  const unsigned short nSpecies = reactiveGas.GetNumSpecies();
  const sunindextype nEq = static_cast<sunindextype>(nSpecies + 1);

  const su2double relTol = min(config->GetChemistry_T_RTL(), config->GetChemistry_Y_RTL());
  const su2double tMin = config->GetChemistry_Clip_TMin();
  const su2double tMax = config->GetChemistry_Clip_TMax();
  const su2double yMin = config->GetChemistry_Clip_YMin();

  struct SubstepDiagnostics {
    bool success = false;
    unsigned long internalSteps = 0;
    int cvodeReturnFlag = 0;
    su2double maxClipCorrection = 0.0;
    su2double sumYCorrection = 0.0;
    su2double maxRenormalizationCorrection = 0.0;
    su2double renormalizationFactor = 1.0;
    su2double temperatureClipCorrection = 0.0;
    std::vector<unsigned short> clippedSpecies;
  };

  auto tryIntegrate = [&](su2double substepDt, su2double& localTemperature,
                          std::vector<su2double>& localMassFractions,
                          SubstepDiagnostics& diagnostics) -> bool {
    SUNContext sunctx = nullptr;
    N_Vector state = nullptr;
    N_Vector absTol = nullptr;
    N_Vector constraints = nullptr;
    SUNMatrix jacobian = nullptr;
    SUNLinearSolver linearSolver = nullptr;
    void* cvodeMem = nullptr;

    auto cleanup = [&]() {
      if (cvodeMem) CVodeFree(&cvodeMem);
      if (linearSolver) SUNLinSolFree(linearSolver);
      if (jacobian) SUNMatDestroy(jacobian);
      if (constraints) N_VDestroy(constraints);
      if (absTol) N_VDestroy(absTol);
      if (state) N_VDestroy(state);
      if (sunctx) SUNContext_Free(&sunctx);
    };

    if (SUNContext_Create(SUN_COMM_NULL, &sunctx) != SUN_SUCCESS) {
      cleanup();
      return false;
    }

    state = N_VNew_Serial(nEq, sunctx);
    absTol = N_VNew_Serial(nEq, sunctx);
    constraints = N_VNew_Serial(nEq, sunctx);
    if (!state || !absTol || !constraints) {
      cleanup();
      return false;
    }

    auto* stateData = N_VGetArrayPointer(state);
    auto* absTolData = N_VGetArrayPointer(absTol);
    auto* constraintData = N_VGetArrayPointer(constraints);

    stateData[0] = localTemperature;
    absTolData[0] = config->GetChemistry_T_ATL();
    constraintData[0] = 1.0;
    for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      stateData[iSpecies + 1] = localMassFractions[iSpecies];
      absTolData[iSpecies + 1] = config->GetChemistry_Y_ATL();
      constraintData[iSpecies + 1] = 1.0;
    }

    ChemistryRHSData rhsData;
    rhsData.reactiveGas = &reactiveGas;
    rhsData.density = density;
    rhsData.massFractions.resize(nSpecies, 0.0);
    rhsData.wdot.resize(nSpecies, 0.0);
    rhsData.speciesEnthalpies.resize(nSpecies, 0.0);
    rhsData.speciesInternalEnergies.resize(nSpecies, 0.0);

    cvodeMem = CVodeCreate(CV_BDF, sunctx);
    if (!cvodeMem) {
      cleanup();
      return false;
    }

    if (!CVodeCallSucceeded(CVodeInit(cvodeMem, DetonationChemistryRHS, 0.0, state)) ||
        !CVodeCallSucceeded(CVodeSVtolerances(cvodeMem, relTol, absTol)) ||
        !CVodeCallSucceeded(CVodeSetUserData(cvodeMem, &rhsData)) ||
        !CVodeCallSucceeded(CVodeSetConstraints(cvodeMem, constraints)) ||
        !CVodeCallSucceeded(CVodeSetMaxNumSteps(cvodeMem, 20000)) ||
        !CVodeCallSucceeded(CVodeSetMaxStep(cvodeMem, substepDt))) {
      cleanup();
      return false;
    }

    jacobian = SUNDenseMatrix(nEq, nEq, sunctx);
    linearSolver = SUNLinSol_Dense(state, jacobian, sunctx);
    if (!jacobian || !linearSolver ||
        !CVodeCallSucceeded(CVodeSetLinearSolver(cvodeMem, linearSolver, jacobian))) {
      cleanup();
      return false;
    }

    sunrealtype tReached = 0.0;
    const int flag = CVode(cvodeMem, substepDt, state, &tReached, CV_NORMAL);
    diagnostics.cvodeReturnFlag = flag;
    if (!CVodeCallSucceeded(flag) || tReached < substepDt) {
      cleanup();
      return false;
    }

    long int numSteps = 0;
    if (CVodeGetNumSteps(cvodeMem, &numSteps) == CV_SUCCESS && numSteps > 0) {
      diagnostics.internalSteps = static_cast<unsigned long>(numSteps);
    }

    const su2double rawTemperature = static_cast<su2double>(stateData[0]);
    diagnostics.temperatureClipCorrection = fabs(min(max(rawTemperature, tMin), tMax) - rawTemperature);
    localTemperature = min(max(rawTemperature, tMin), tMax);

    std::vector<su2double> clippedMassFractions(nSpecies, 0.0);
    su2double sumClipped = 0.0;
    su2double sumRaw = 0.0;
    for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      const su2double rawValue = static_cast<su2double>(stateData[iSpecies + 1]);
      sumRaw += rawValue;
      const su2double clippedValue = max(rawValue, yMin);
      clippedMassFractions[iSpecies] = clippedValue;
      sumClipped += clippedValue;
      diagnostics.maxClipCorrection = max(diagnostics.maxClipCorrection, fabs(clippedValue - rawValue));
      if (fabs(clippedValue - rawValue) > 0.0) diagnostics.clippedSpecies.push_back(iSpecies);
    }
    diagnostics.sumYCorrection = fabs(sumRaw - 1.0);
    diagnostics.renormalizationFactor = sumClipped;

    if (sumClipped <= 0.0) {
      localMassFractions = clippedMassFractions;
      reactiveGas.NormalizeMassFractions(localMassFractions, yMin);
    } else {
      for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
        const su2double normalizedValue = clippedMassFractions[iSpecies] / sumClipped;
        diagnostics.maxRenormalizationCorrection =
            max(diagnostics.maxRenormalizationCorrection, fabs(normalizedValue - clippedMassFractions[iSpecies]));
        localMassFractions[iSpecies] = normalizedValue;
      }
    }

    diagnostics.success = true;
    cleanup();
    return true;
  };

  const unsigned long maxRetrySubsteps = max(config->GetChemistry_Max_Substeps(), 1ul);
  unsigned long nRetrySubsteps = 1;
  for (; nRetrySubsteps <= maxRetrySubsteps; nRetrySubsteps *= 2) {
    su2double candidateTemperature = result.temperature;
    std::vector<su2double> candidateMassFractions = result.massFractions;
    unsigned long internalStepsTotal = 0;
    int lastCVodeFlag = 0;
    su2double maxClipCorrection = 0.0;
    su2double maxSumYCorrection = 0.0;
    su2double maxRenormCorrection = 0.0;
    su2double maxTemperatureClipCorrection = 0.0;
    su2double finalRenormalizationFactor = 1.0;
    std::vector<bool> clippedMask(nSpecies, false);

    bool success = true;
    const su2double substepDt = dt / static_cast<su2double>(nRetrySubsteps);
    for (unsigned long iSubstep = 0; iSubstep < nRetrySubsteps; ++iSubstep) {
      SubstepDiagnostics diagnostics;
      success = tryIntegrate(substepDt, candidateTemperature, candidateMassFractions, diagnostics);
      internalStepsTotal += diagnostics.internalSteps;
      lastCVodeFlag = diagnostics.cvodeReturnFlag;
      maxClipCorrection = max(maxClipCorrection, diagnostics.maxClipCorrection);
      maxSumYCorrection = max(maxSumYCorrection, diagnostics.sumYCorrection);
      maxRenormCorrection = max(maxRenormCorrection, diagnostics.maxRenormalizationCorrection);
      maxTemperatureClipCorrection = max(maxTemperatureClipCorrection, diagnostics.temperatureClipCorrection);
      finalRenormalizationFactor = diagnostics.renormalizationFactor;
      for (auto iSpecies : diagnostics.clippedSpecies) clippedMask[iSpecies] = true;
      if (!success) break;
    }

    if (!success) continue;

    result.success = true;
    result.substeps = nRetrySubsteps;
    result.internalSteps = internalStepsTotal;
    result.cvodeReturnFlag = lastCVodeFlag;
    result.temperature = candidateTemperature;
    result.massFractions = std::move(candidateMassFractions);
    result.maxClipCorrection = maxClipCorrection;
    result.sumYCorrection = maxSumYCorrection;
    result.maxRenormalizationCorrection = maxRenormCorrection;
    result.renormalizationFactor = finalRenormalizationFactor;
    result.temperatureClipCorrection = maxTemperatureClipCorrection;
    result.clippedSpecies.clear();
    for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      if (clippedMask[iSpecies]) result.clippedSpecies.push_back(iSpecies);
    }

    std::vector<su2double> wdot(nSpecies, 0.0), speciesEnthalpies(nSpecies, 0.0), speciesInternalEnergies(nSpecies, 0.0);
    su2double cvMix = 0.0;
    reactiveGas.ComputeChemistrySourceTerms(density, result.temperature, result.massFractions, wdot, speciesEnthalpies,
                                            speciesInternalEnergies, cvMix, result.qdot);
    result.maxAbsWdot = 0.0;
    for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      result.maxAbsWdot = max(result.maxAbsWdot, fabs(wdot[iSpecies]));
    }
    for (unsigned short iRank = 0; iRank < result.dominantSpecies.size(); ++iRank) {
      su2double bestMagnitude = -1.0;
      int bestSpecies = -1;
      su2double bestWdot = 0.0;
      for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
        bool alreadySelected = false;
        for (unsigned short j = 0; j < iRank; ++j) {
          if (result.dominantSpecies[j] == static_cast<int>(iSpecies)) {
            alreadySelected = true;
            break;
          }
        }
        if (alreadySelected) continue;
        const su2double magnitude = fabs(wdot[iSpecies]);
        if (magnitude > bestMagnitude) {
          bestMagnitude = magnitude;
          bestSpecies = static_cast<int>(iSpecies);
          bestWdot = wdot[iSpecies];
        }
      }
      result.dominantSpecies[iRank] = bestSpecies;
      result.dominantWdot[iRank] = bestWdot;
    }

    result.sumYAfter = 0.0;
    result.minYAfter = std::numeric_limits<su2double>::max();
    for (auto value : result.massFractions) {
      result.sumYAfter += value;
      result.minYAfter = min(result.minYAfter, value);
    }
    if (result.minYAfter == std::numeric_limits<su2double>::max()) result.minYAfter = 0.0;

    evaluateState(result.temperature, result.massFractions, result.internalEnergyAfter, result.pressureAfter);
    result.chemistryEnergyAbsError = result.internalEnergyAfter - eOld;
    result.chemistryEnergyRelError =
        fabs(result.chemistryEnergyAbsError) / max(fabs(eOld), su2double(1.0e-30));
    return result;
  }

  result.success = false;
  result.substeps = maxRetrySubsteps;
  result.cvodeReturnFlag = CV_TOO_MUCH_WORK;
  result.qdot = 0.0;
  return result;
#endif
}
