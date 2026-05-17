#pragma once

#include "../CNumerics.hpp"

class CDetonationNumerics : public CNumerics {
protected:
  unsigned short nSpecies = 0;
  unsigned short nPrimVar = 0;
  unsigned short RHO_INDEX = 0;
  unsigned short T_INDEX = 0;
  unsigned short VEL_INDEX = 0;
  unsigned short P_INDEX = 0;
  unsigned short H_INDEX = 0;
  unsigned short A_INDEX = 0;
  su2double* Flux = nullptr;

  void GetInviscidProjFlux(const su2double* conserved,
                           const su2double* primitive,
                           const su2double* normal,
                           su2double* projectedFlux) const;

public:
  CDetonationNumerics(unsigned short nDim, unsigned short nVar, unsigned short nPrimVar, const CConfig* config);
  ~CDetonationNumerics() override;
};

class CUpwLLF_Detonation final : public CDetonationNumerics {
public:
  CUpwLLF_Detonation(unsigned short nDim, unsigned short nVar, unsigned short nPrimVar, const CConfig* config);
  ResidualType<> ComputeResidual(const CConfig* config) override;
};

class CUpwHLLC_Detonation final : public CDetonationNumerics {
private:
  su2double* FluxL = nullptr;
  su2double* FluxR = nullptr;
  su2double* UStar = nullptr;

  static unsigned long FallbackCount;
  static unsigned long NegativeStateCount;
  static su2double MinPressure;
  static su2double MinDensity;

  void ComputeLLFFallback();
  bool ComputeStarState(const su2double* conserved, const su2double* primitive,
                        const su2double* unitNormal, su2double waveSpeed,
                        su2double contactSpeed, su2double normalVelocity,
                        su2double* starState, su2double& starPressure,
                        su2double& starDensity) const;

public:
  CUpwHLLC_Detonation(unsigned short nDim, unsigned short nVar, unsigned short nPrimVar, const CConfig* config);
  ~CUpwHLLC_Detonation() override;

  ResidualType<> ComputeResidual(const CConfig* config) override;

  static void ResetDiagnostics();
  static unsigned long GetFallbackCount();
  static unsigned long GetNegativeStateCount();
  static su2double GetMinPressure();
  static su2double GetMinDensity();
};
