#include "../../include/fluid/CCanteraReactiveGas.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <utility>

#ifdef HAVE_CANTERA
#include "cantera/base/Solution.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/thermo/ThermoPhase.h"
#endif

CCanteraReactiveGas::CCanteraReactiveGas(const CConfig* config, unsigned short nDimIn)
    : nDim(nDimIn),
      nSpecies(config->GetnSpecies()),
      MechanismFile(config->GetChemistry_File()) {
  CurrentMassFractions.assign(nSpecies, 0.0);
  CurrentSpeciesDensities.assign(nSpecies, 0.0);
  MolecularWeights.assign(nSpecies, 0.0);
  SpeciesEnthalpies.assign(nSpecies, 0.0);
  SpeciesInternalEnergies.assign(nSpecies, 0.0);

  if (nSpecies == 0) {
    SU2_MPI::Error("CANTERA_REACTIVE requires at least one species.", CURRENT_FUNCTION);
  }

#ifndef HAVE_CANTERA
  SU2_MPI::Error("DETONATION_EULER was selected but this binary was not compiled with Cantera support. Reconfigure with -Denable-cantera=true.", CURRENT_FUNCTION);
#else
  try {
    CanteraSolution = Cantera::newSolution(MechanismFile, "", "none");
    Thermo = CanteraSolution->thermo();
    KineticsModel = CanteraSolution->kinetics();
  } catch (const std::exception& err) {
    SU2_MPI::Error("Failed to load Cantera mechanism '" + MechanismFile + "': " + err.what(), CURRENT_FUNCTION);
  }

  if (!Thermo || !KineticsModel) {
    SU2_MPI::Error("Cantera mechanism initialization did not provide thermo/kinetics objects.", CURRENT_FUNCTION);
  }

  if (Thermo->nSpecies() != nSpecies) {
    std::ostringstream msg;
    msg << "Phase-1 DETONATION_EULER currently requires GAS_COMPOSITION to provide one entry for each species in the Cantera mechanism."
        << " Mechanism '" << MechanismFile << "' contains " << Thermo->nSpecies()
        << " species, but GAS_COMPOSITION provided " << nSpecies << ".";
    SU2_MPI::Error(msg.str(), CURRENT_FUNCTION);
  }

  for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    MolecularWeights[iSpecies] = Thermo->molecularWeight(iSpecies);
  }

  SetMassFractions(config->GetGas_Composition(), nSpecies);
  const su2double pressure0 = (config->GetPressure_FreeStreamND() > 0.0)
                                  ? config->GetPressure_FreeStreamND()
                                  : config->GetPressure_FreeStream();
  const su2double temperature0 = (config->GetTemperature_FreeStreamND() > 0.0)
                                     ? config->GetTemperature_FreeStreamND()
                                     : config->GetTemperature_FreeStream();
  SetTDState_PT(pressure0, temperature0);

  if (SU2_MPI::GetRank() == MASTER_NODE) {
    cout << "Detonation chemistry backend: Cantera mechanism '" << MechanismFile
         << "' loaded with " << nSpecies << " species." << endl;
  }
#endif
}

void CCanteraReactiveGas::NormalizeMassFractions(std::vector<su2double>& massFractions, su2double ymin) const {
  su2double sum = 0.0;
  for (auto& value : massFractions) {
    value = max(value, ymin);
    sum += value;
  }
  if (sum <= 0.0) {
    const su2double uniform = 1.0 / static_cast<su2double>(massFractions.size());
    for (auto& value : massFractions) value = uniform;
    return;
  }
  for (auto& value : massFractions) value /= sum;
}

void CCanteraReactiveGas::SetMassFractions(const su2double* massFractions, unsigned short count) {
  std::fill(CurrentMassFractions.begin(), CurrentMassFractions.end(), 0.0);
  for (unsigned short iSpecies = 0; iSpecies < min(count, nSpecies); ++iSpecies) {
    CurrentMassFractions[iSpecies] = massFractions[iSpecies];
  }
  NormalizeMassFractions(CurrentMassFractions, 0.0);
}

void CCanteraReactiveGas::SetSpeciesDensities(const su2double* speciesDensities, unsigned short count) {
  Density = 0.0;
  std::fill(CurrentSpeciesDensities.begin(), CurrentSpeciesDensities.end(), 0.0);
  for (unsigned short iSpecies = 0; iSpecies < min(count, nSpecies); ++iSpecies) {
    CurrentSpeciesDensities[iSpecies] = max(speciesDensities[iSpecies], 0.0);
    Density += CurrentSpeciesDensities[iSpecies];
  }

  if (Density <= 0.0) Density = 1.0e-16;

  for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    CurrentMassFractions[iSpecies] = CurrentSpeciesDensities[iSpecies] / Density;
  }
  NormalizeMassFractions(CurrentMassFractions, 0.0);
}

void CCanteraReactiveGas::ApplyMassFractionsToBackend() {
#ifdef HAVE_CANTERA
  Thermo->setMassFractions_NoNorm(CurrentMassFractions.data());
#endif
}

void CCanteraReactiveGas::UpdateStateFromBackend() {
#ifdef HAVE_CANTERA
  Thermo->getMassFractions(CurrentMassFractions.data());

  Density = Thermo->density();
  Pressure = Thermo->pressure();
  Temperature = Thermo->temperature();
  StaticEnergy = Thermo->intEnergy_mass();
  Entropy = Thermo->entropy_mass();
  Cp = Thermo->cp_mass();
  Cv = Thermo->cv_mass();

  const su2double meanMolWeight = Thermo->meanMolecularWeight();
  MixtureGasConstant = Cantera::GasConstant / max(meanMolWeight, su2double(1.0e-16));
  MixtureGamma = (Cv > 0.0) ? Cp / Cv : 1.4;

  su2double soundSpeed = 0.0;
  try {
    soundSpeed = Thermo->soundSpeed();
  } catch (const std::exception&) {
    soundSpeed = 0.0;
  }

  if (!(soundSpeed > 0.0) && Density > 0.0) {
    soundSpeed = sqrt(max(MixtureGamma * Pressure / Density, su2double(0.0)));
  }
  SoundSpeed2 = soundSpeed * soundSpeed;

  if (Density > 0.0) {
    dPdrho_e = Pressure / Density;
    dPde_rho = Density * MixtureGasConstant / max(Cv, su2double(1.0e-16));
  } else {
    dPdrho_e = 0.0;
    dPde_rho = 0.0;
  }
  dTdrho_e = 0.0;
  dTde_rho = 1.0 / max(Cv, su2double(1.0e-16));

  for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    CurrentSpeciesDensities[iSpecies] = Density * CurrentMassFractions[iSpecies];
  }
#endif
}

void CCanteraReactiveGas::SetCanteraState_rhoT(su2double rho, su2double temperature) {
#ifdef HAVE_CANTERA
  ApplyMassFractionsToBackend();
  Thermo->setState_TD(temperature, rho);
  UpdateStateFromBackend();
#else
  Density = rho;
  Temperature = temperature;
#endif
}

void CCanteraReactiveGas::SetTDState_rhoe(su2double rho, su2double e) {
#ifdef HAVE_CANTERA
  try {
    ApplyMassFractionsToBackend();
    Thermo->setState_UV(e, 1.0 / max(rho, su2double(1.0e-16)));
    UpdateStateFromBackend();
  } catch (const std::exception& err) {
    SU2_MPI::Error("Cantera failed to recover state from (rho,e): " + std::string(err.what()), CURRENT_FUNCTION);
  }
#else
  Density = rho;
  StaticEnergy = e;
#endif
}

void CCanteraReactiveGas::SetTDState_PT(su2double pressure, su2double temperature) {
#ifdef HAVE_CANTERA
  try {
    Thermo->setState_TPY(temperature, pressure, CurrentMassFractions.data());
    UpdateStateFromBackend();
  } catch (const std::exception& err) {
    SU2_MPI::Error("Cantera failed to set state from (P,T): " + std::string(err.what()), CURRENT_FUNCTION);
  }
#else
  Pressure = pressure;
  Temperature = temperature;
#endif
}

void CCanteraReactiveGas::SetTDState_Prho(su2double pressure, su2double rho) {
#ifdef HAVE_CANTERA
  try {
    ApplyMassFractionsToBackend();
    const su2double meanMolWeight = Thermo->meanMolecularWeight();
    const su2double gasConstant = Cantera::GasConstant / max(meanMolWeight, su2double(1.0e-16));
    const su2double temperature = pressure / max(rho * gasConstant, su2double(1.0e-16));
    Thermo->setState_TD(temperature, rho);
    UpdateStateFromBackend();
  } catch (const std::exception& err) {
    SU2_MPI::Error("Cantera failed to set state from (P,rho): " + std::string(err.what()), CURRENT_FUNCTION);
  }
#else
  Pressure = pressure;
  Density = rho;
#endif
}

void CCanteraReactiveGas::SetTDState_rhoT(su2double rho, su2double temperature) {
#ifdef HAVE_CANTERA
  try {
    SetCanteraState_rhoT(rho, temperature);
  } catch (const std::exception& err) {
    SU2_MPI::Error("Cantera failed to set state from (rho,T): " + std::string(err.what()), CURRENT_FUNCTION);
  }
#else
  Density = rho;
  Temperature = temperature;
#endif
}

void CCanteraReactiveGas::ComputeChemistryThermo(su2double temperature,
                                                 const std::vector<su2double>& massFractions,
                                                 std::vector<su2double>& speciesEnthalpiesOut,
                                                 std::vector<su2double>& speciesInternalEnergiesOut,
                                                 su2double& cvMix) {
#ifdef HAVE_CANTERA
  std::vector<su2double> localMassFractions = massFractions;
  NormalizeMassFractions(localMassFractions, 0.0);
  SetMassFractions(localMassFractions.data(), nSpecies);
  SetCanteraState_rhoT(1.0, temperature);

  Thermo->getPartialMolarEnthalpies(SpeciesEnthalpies.data());
  Thermo->getPartialMolarIntEnergies(SpeciesInternalEnergies.data());

  speciesEnthalpiesOut.resize(nSpecies);
  speciesInternalEnergiesOut.resize(nSpecies);
  for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    const su2double molecularWeight = max(MolecularWeights[iSpecies], su2double(1.0e-30));
    speciesEnthalpiesOut[iSpecies] = SpeciesEnthalpies[iSpecies] / molecularWeight;
    speciesInternalEnergiesOut[iSpecies] = SpeciesInternalEnergies[iSpecies] / molecularWeight;
  }
  cvMix = Thermo->cv_mass();
#else
  cvMix = MixtureGasConstant / max(MixtureGamma - 1.0, su2double(1.0e-16));
  speciesEnthalpiesOut.assign(nSpecies, (cvMix + MixtureGasConstant) * temperature);
  speciesInternalEnergiesOut.assign(nSpecies, cvMix * temperature);
#endif
}

void CCanteraReactiveGas::ComputeChemistrySourceTerms(su2double density,
                                                      su2double temperature,
                                                      const std::vector<su2double>& massFractions,
                                                      std::vector<su2double>& wdot,
                                                      std::vector<su2double>& speciesEnthalpiesOut,
                                                      std::vector<su2double>& speciesInternalEnergiesOut,
                                                      su2double& cvMix,
                                                      su2double& qdot) {
#ifdef HAVE_CANTERA
  std::vector<su2double> localMassFractions = massFractions;
  NormalizeMassFractions(localMassFractions, 0.0);
  SetMassFractions(localMassFractions.data(), nSpecies);
  SetCanteraState_rhoT(density, temperature);

  Thermo->getPartialMolarEnthalpies(SpeciesEnthalpies.data());
  Thermo->getPartialMolarIntEnergies(SpeciesInternalEnergies.data());

  wdot.resize(nSpecies);
  speciesEnthalpiesOut.resize(nSpecies);
  speciesInternalEnergiesOut.resize(nSpecies);
  std::vector<double> wdotMolar(nSpecies, 0.0);
  KineticsModel->getNetProductionRates(wdotMolar.data());

  cvMix = Thermo->cv_mass();
  qdot = 0.0;
  for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    const su2double molecularWeight = max(MolecularWeights[iSpecies], su2double(1.0e-30));
    wdot[iSpecies] = wdotMolar[iSpecies] * molecularWeight;
    speciesEnthalpiesOut[iSpecies] = SpeciesEnthalpies[iSpecies] / molecularWeight;
    speciesInternalEnergiesOut[iSpecies] = SpeciesInternalEnergies[iSpecies] / molecularWeight;
    qdot -= speciesEnthalpiesOut[iSpecies] * wdot[iSpecies];
  }
#else
  ComputeChemistryThermo(temperature, massFractions, speciesEnthalpiesOut, speciesInternalEnergiesOut, cvMix);
  wdot.assign(nSpecies, 0.0);
  qdot = 0.0;
#endif
}
