#include "../../include/variables/CDetonationEulerVariable.hpp"

CDetonationEulerVariable::CDetonationEulerVariable(su2double pressure,
                                                   const su2double* massFractionsIn,
                                                   const su2double* velocity,
                                                   su2double temperature,
                                                   unsigned long nPointIn,
                                                   unsigned long nDimIn,
                                                   unsigned long nVarIn,
                                                   unsigned long nPrimVarIn,
                                                   unsigned long nPrimVarGradIn,
                                                   const CConfig* config,
                                                   CCanteraReactiveGas* fluidModel)
    : CFlowVariable(nPointIn, nDimIn, nVarIn, nPrimVarIn, nPrimVarGradIn, config),
      indices(nDimIn, config->GetnSpecies()) {
  nSpecies = config->GetnSpecies();
  nSecondaryVar = 2;
  nSecondaryVarGrad = 2;

  Secondary.resize(nPoint, nSecondaryVar) = su2double(0.0);
  MassFractions.resize(nPoint, nSpecies) = su2double(0.0);
  MixtureGamma.resize(nPoint) = su2double(0.0);

  fluidModel->SetMassFractions(massFractionsIn, nSpecies);
  fluidModel->SetTDState_PT(pressure, temperature);

  const su2double rho = fluidModel->GetDensity();
  su2double velocity2 = 0.0;
  for (unsigned short iDim = 0; iDim < nDim; ++iDim) velocity2 += velocity[iDim] * velocity[iDim];
  const su2double energy = fluidModel->GetStaticEnergy() + 0.5 * velocity2;

  for (unsigned long iPoint = 0; iPoint < nPoint; ++iPoint) {
    for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      Solution(iPoint, iSpecies) = rho * massFractionsIn[iSpecies];
      MassFractions(iPoint, iSpecies) = massFractionsIn[iSpecies];
      Primitive(iPoint, iSpecies) = rho * massFractionsIn[iSpecies];
    }
    for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
      Solution(iPoint, nSpecies + iDim) = rho * velocity[iDim];
      Primitive(iPoint, indices.Velocity() + iDim) = velocity[iDim];
    }
    Solution(iPoint, nVar - 1) = rho * energy;
  }

  Solution_Old = Solution;
}

bool CDetonationEulerVariable::SetPrimVar(unsigned long iPoint, CFluidModel* fluidModelIn) {
  auto* reactiveGas = static_cast<CCanteraReactiveGas*>(fluidModelIn);

  bool physical = true;
  if (SetDensity(iPoint)) physical = false;

  SetVelocity(iPoint);

  std::vector<su2double> speciesDensities(nSpecies, 0.0);
  for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    speciesDensities[iSpecies] = max(Solution(iPoint, iSpecies), 0.0);
    Primitive(iPoint, iSpecies) = speciesDensities[iSpecies];
  }

  reactiveGas->SetSpeciesDensities(speciesDensities.data(), nSpecies);

  const su2double rho = Primitive(iPoint, indices.Density());
  const su2double totalEnergy = Solution(iPoint, nVar - 1) / rho;
  const su2double staticEnergy = totalEnergy - 0.5 * Velocity2(iPoint);
  reactiveGas->SetTDState_rhoe(rho, staticEnergy);

  physical = physical && !SetPressure(iPoint, reactiveGas->GetPressure());
  physical = physical && !SetSoundSpeed(iPoint, reactiveGas->GetSoundSpeed2());
  physical = physical && !SetTemperature(iPoint, reactiveGas->GetTemperature());
  SetEnthalpy(iPoint);

  const auto& normalizedMassFractions = reactiveGas->GetMassFractions();
  for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    MassFractions(iPoint, iSpecies) = normalizedMassFractions[iSpecies];
  }
  MixtureGamma(iPoint) = reactiveGas->GetGamma();

  if (!physical) {
    for (unsigned long iVar = 0; iVar < nVar; ++iVar) Solution(iPoint, iVar) = Solution_Old(iPoint, iVar);
  }

  return physical;
}

void CDetonationEulerVariable::SetSecondaryVar(unsigned long iPoint, CFluidModel* fluidModel) {
  SetdPdrho_e(iPoint, fluidModel->GetdPdrho_e());
  SetdPde_rho(iPoint, fluidModel->GetdPde_rho());
}
