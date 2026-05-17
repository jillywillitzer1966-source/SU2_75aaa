#pragma once

#include <memory>
#include <string>
#include <vector>

#include "CFluidModel.hpp"

namespace Cantera {
class Solution;
class ThermoPhase;
class Kinetics;
}

class CCanteraReactiveGas final : public CFluidModel {
private:
  unsigned short nDim = 0;
  unsigned short nSpecies = 0;
  su2double MixtureGamma = 1.4;
  su2double MixtureGasConstant = 287.058;
  std::vector<su2double> CurrentMassFractions;
  std::vector<su2double> CurrentSpeciesDensities;
  std::vector<su2double> MolecularWeights;
  std::vector<su2double> SpeciesEnthalpies;
  std::vector<su2double> SpeciesInternalEnergies;
  std::string MechanismFile;

#ifdef HAVE_CANTERA
  std::shared_ptr<Cantera::Solution> CanteraSolution;
  std::shared_ptr<Cantera::ThermoPhase> Thermo;
  std::shared_ptr<Cantera::Kinetics> KineticsModel;
#endif

  void UpdateStateFromBackend();
  void ApplyMassFractionsToBackend();
  void SetCanteraState_rhoT(su2double rho, su2double temperature);

public:
  explicit CCanteraReactiveGas(const CConfig* config, unsigned short nDimIn);

  unsigned short GetNumSpecies() const { return nSpecies; }

  const std::vector<su2double>& GetMassFractions() const { return CurrentMassFractions; }

  void NormalizeMassFractions(std::vector<su2double>& massFractions, su2double ymin) const;

  void SetMassFractions(const su2double* massFractions, unsigned short count);

  void SetSpeciesDensities(const su2double* speciesDensities, unsigned short count);

  void SetTDState_rhoe(su2double rho, su2double e) override;

  void SetTDState_PT(su2double pressure, su2double temperature) override;

  void SetTDState_Prho(su2double pressure, su2double rho) override;

  void SetTDState_rhoT(su2double rho, su2double temperature) override;

  void ComputeChemistryThermo(su2double temperature,
                              const std::vector<su2double>& massFractions,
                              std::vector<su2double>& speciesEnthalpies,
                              std::vector<su2double>& speciesInternalEnergies,
                              su2double& cvMix);

  void ComputeChemistrySourceTerms(su2double density,
                                   su2double temperature,
                                   const std::vector<su2double>& massFractions,
                                   std::vector<su2double>& wdot,
                                   std::vector<su2double>& speciesEnthalpies,
                                   std::vector<su2double>& speciesInternalEnergies,
                                   su2double& cvMix,
                                   su2double& qdot);

  su2double GetGamma() const { return MixtureGamma; }
};
