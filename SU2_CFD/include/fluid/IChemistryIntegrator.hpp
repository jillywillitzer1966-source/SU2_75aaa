#pragma once

#include <array>
#include <vector>

#include "../../../Common/include/basic_types/datatype_structure.hpp"

struct ChemistryAdvanceResult {
  bool success = false;
  unsigned long substeps = 0;
  unsigned long internalSteps = 0;
  int cvodeReturnFlag = 0;
  su2double temperature = 0.0;
  su2double qdot = 0.0;
  su2double chemistryEnergyAbsError = 0.0;
  su2double chemistryEnergyRelError = 0.0;
  su2double density = 0.0;
  su2double temperatureBefore = 0.0;
  su2double pressureBefore = 0.0;
  su2double pressureAfter = 0.0;
  su2double internalEnergyBefore = 0.0;
  su2double internalEnergyAfter = 0.0;
  su2double sumYBefore = 0.0;
  su2double sumYAfter = 0.0;
  su2double minYBefore = 0.0;
  su2double minYAfter = 0.0;
  su2double maxClipCorrection = 0.0;
  su2double sumYCorrection = 0.0;
  su2double maxRenormalizationCorrection = 0.0;
  su2double renormalizationFactor = 1.0;
  su2double temperatureClipCorrection = 0.0;
  su2double maxAbsWdot = 0.0;
  std::array<int, 3> dominantSpecies{{-1, -1, -1}};
  std::array<su2double, 3> dominantWdot{{0.0, 0.0, 0.0}};
  std::vector<unsigned short> clippedSpecies;
  std::vector<su2double> massFractions;
};

class IChemistryIntegrator {
public:
  virtual ~IChemistryIntegrator() = default;

  virtual ChemistryAdvanceResult IntegrateCell(su2double dt,
                                               su2double density,
                                               su2double temperature,
                                               const std::vector<su2double>& massFractions) = 0;
};
