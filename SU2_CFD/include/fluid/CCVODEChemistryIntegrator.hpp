#pragma once

#include "IChemistryIntegrator.hpp"

class CConfig;
class CCanteraReactiveGas;

class CCVODEChemistryIntegrator final : public IChemistryIntegrator {
private:
  CCanteraReactiveGas& reactiveGas;
  const CConfig* config = nullptr;

public:
  CCVODEChemistryIntegrator(CCanteraReactiveGas& reactiveGasIn, const CConfig* configIn);

  ChemistryAdvanceResult IntegrateCell(su2double dt,
                                       su2double density,
                                       su2double temperature,
                                       const std::vector<su2double>& massFractions) override;
};
