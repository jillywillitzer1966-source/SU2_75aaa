#include "../../../include/numerics/detonation/CDetonationNumerics.hpp"

#include "../../../../Common/include/toolboxes/geometry_toolbox.hpp"

#include <algorithm>
#include <limits>

CDetonationNumerics::CDetonationNumerics(unsigned short nDimIn, unsigned short nVarIn, unsigned short nPrimVarIn,
                                         const CConfig* config)
    : CNumerics(nDimIn, nVarIn, config) {
  nSpecies = nVar - nDim - 1;
  nPrimVar = nPrimVarIn;
  RHO_INDEX = nSpecies + nDim + 2;
  T_INDEX = nSpecies;
  VEL_INDEX = nSpecies + 1;
  P_INDEX = nSpecies + nDim + 1;
  H_INDEX = nSpecies + nDim + 3;
  A_INDEX = nSpecies + nDim + 4;
  Flux = new su2double[nVar]();
}

CDetonationNumerics::~CDetonationNumerics() { delete[] Flux; }

void CDetonationNumerics::GetInviscidProjFlux(const su2double* conserved,
                                              const su2double* primitive,
                                              const su2double* normal,
                                              su2double* projectedFlux) const {
  for (unsigned short iVar = 0; iVar < nVar; ++iVar) projectedFlux[iVar] = 0.0;

  const su2double pressure = primitive[P_INDEX];
  const su2double enthalpy = primitive[H_INDEX];
  su2double projVelocity = 0.0;
  for (unsigned short iDim = 0; iDim < nDim; ++iDim) projVelocity += primitive[VEL_INDEX + iDim] * normal[iDim];

  for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    projectedFlux[iSpecies] = conserved[iSpecies] * projVelocity;
  }

  for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
    projectedFlux[nSpecies + iDim] = conserved[nSpecies + iDim] * projVelocity + pressure * normal[iDim];
  }

  projectedFlux[nVar - 1] = primitive[RHO_INDEX] * projVelocity * enthalpy;
}

CUpwLLF_Detonation::CUpwLLF_Detonation(unsigned short nDim, unsigned short nVar,
                                       unsigned short nPrimVar, const CConfig* config)
    : CDetonationNumerics(nDim, nVar, nPrimVar, config) {
  if (SU2_MPI::GetRank() == MASTER_NODE) {
    cout << "DETONATION_EULER convective flux: DETONATION_LLF selected. "
         << "This path uses a local Lax-Friedrichs/Rusanov-like placeholder flux; true AUSM+UP2 is not implemented."
         << endl;
  }
}

CNumerics::ResidualType<> CUpwLLF_Detonation::ComputeResidual(const CConfig*) {
  su2double fluxI[256] = {0.0};
  su2double fluxJ[256] = {0.0};

  GetInviscidProjFlux(U_i, V_i, Normal, fluxI);
  GetInviscidProjFlux(U_j, V_j, Normal, fluxJ);

  const su2double projVelocityI = GeometryToolbox::DotProduct(nDim, &V_i[VEL_INDEX], Normal);
  const su2double projVelocityJ = GeometryToolbox::DotProduct(nDim, &V_j[VEL_INDEX], Normal);
  const su2double lambda = max(fabs(projVelocityI) + V_i[A_INDEX] * GeometryToolbox::Norm(nDim, Normal),
                               fabs(projVelocityJ) + V_j[A_INDEX] * GeometryToolbox::Norm(nDim, Normal));

  for (unsigned short iVar = 0; iVar < nVar; ++iVar) {
    Flux[iVar] = 0.5 * (fluxI[iVar] + fluxJ[iVar]) - 0.5 * lambda * (U_j[iVar] - U_i[iVar]);
  }

  return ResidualType<>(Flux, nullptr, nullptr);
}

unsigned long CUpwHLLC_Detonation::FallbackCount = 0;
unsigned long CUpwHLLC_Detonation::NegativeStateCount = 0;
su2double CUpwHLLC_Detonation::MinPressure = std::numeric_limits<su2double>::max();
su2double CUpwHLLC_Detonation::MinDensity = std::numeric_limits<su2double>::max();

CUpwHLLC_Detonation::CUpwHLLC_Detonation(unsigned short nDim, unsigned short nVar,
                                         unsigned short nPrimVar, const CConfig* config)
    : CDetonationNumerics(nDim, nVar, nPrimVar, config) {
  FluxL = new su2double[nVar]();
  FluxR = new su2double[nVar]();
  UStar = new su2double[nVar]();

  if (SU2_MPI::GetRank() == MASTER_NODE) {
    cout << "DETONATION_EULER convective flux: DETONATION_HLLC selected. "
         << "This candidate HLLC flux keeps chemistry, splitting, and state layout unchanged." << endl;
  }
}

CUpwHLLC_Detonation::~CUpwHLLC_Detonation() {
  delete[] FluxL;
  delete[] FluxR;
  delete[] UStar;
}

void CUpwHLLC_Detonation::ResetDiagnostics() {
  FallbackCount = 0;
  NegativeStateCount = 0;
  MinPressure = std::numeric_limits<su2double>::max();
  MinDensity = std::numeric_limits<su2double>::max();
}

unsigned long CUpwHLLC_Detonation::GetFallbackCount() { return FallbackCount; }
unsigned long CUpwHLLC_Detonation::GetNegativeStateCount() { return NegativeStateCount; }
su2double CUpwHLLC_Detonation::GetMinPressure() {
  return (MinPressure == std::numeric_limits<su2double>::max()) ? 0.0 : MinPressure;
}
su2double CUpwHLLC_Detonation::GetMinDensity() {
  return (MinDensity == std::numeric_limits<su2double>::max()) ? 0.0 : MinDensity;
}

void CUpwHLLC_Detonation::ComputeLLFFallback() {
  su2double fluxI[256] = {0.0};
  su2double fluxJ[256] = {0.0};

  GetInviscidProjFlux(U_i, V_i, Normal, fluxI);
  GetInviscidProjFlux(U_j, V_j, Normal, fluxJ);

  const su2double projVelocityI = GeometryToolbox::DotProduct(nDim, &V_i[VEL_INDEX], Normal);
  const su2double projVelocityJ = GeometryToolbox::DotProduct(nDim, &V_j[VEL_INDEX], Normal);
  const su2double lambda = std::max(fabs(projVelocityI) + V_i[A_INDEX] * GeometryToolbox::Norm(nDim, Normal),
                                    fabs(projVelocityJ) + V_j[A_INDEX] * GeometryToolbox::Norm(nDim, Normal));

  for (unsigned short iVar = 0; iVar < nVar; ++iVar) {
    Flux[iVar] = 0.5 * (fluxI[iVar] + fluxJ[iVar]) - 0.5 * lambda * (U_j[iVar] - U_i[iVar]);
  }
  FallbackCount++;
}

bool CUpwHLLC_Detonation::ComputeStarState(const su2double* conserved, const su2double* primitive,
                                           const su2double* unitNormal, su2double waveSpeed,
                                           su2double contactSpeed, su2double normalVelocity,
                                           su2double* starState, su2double& starPressure,
                                           su2double& starDensity) const {
  const su2double rho = primitive[RHO_INDEX];
  const su2double pressure = primitive[P_INDEX];
  const su2double denom = waveSpeed - contactSpeed;
  const su2double acousticJump = waveSpeed - normalVelocity;

  if (rho <= 0.0 || fabs(denom) < 1.0e-14 || fabs(acousticJump) < 1.0e-14) return false;

  const su2double factor = acousticJump / denom;
  starDensity = rho * factor;

  for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    const su2double massFraction = conserved[iSpecies] / rho;
    starState[iSpecies] = starDensity * massFraction;
  }

  su2double starVelocity2 = 0.0;
  for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
    const su2double velocityStar = primitive[VEL_INDEX + iDim] + (contactSpeed - normalVelocity) * unitNormal[iDim];
    starState[nSpecies + iDim] = starDensity * velocityStar;
    starVelocity2 += velocityStar * velocityStar;
  }

  const su2double specificEnergy = conserved[nVar - 1] / rho;
  const su2double specificStarEnergy =
      specificEnergy + (contactSpeed - normalVelocity) *
                       (contactSpeed + pressure / (rho * acousticJump));
  starState[nVar - 1] = starDensity * specificStarEnergy;

  starPressure = pressure + rho * acousticJump * (contactSpeed - normalVelocity);

  MinPressure = std::min(MinPressure, starPressure);
  MinDensity = std::min(MinDensity, starDensity);

  if (starDensity <= 0.0 || starPressure <= 0.0) return false;
  /*
   * Do not reject states solely because rhoE is below kinetic energy.
   * CCanteraReactiveGas uses Cantera thermochemistry with formation-energy
   * references, so the mixture internal energy may be negative while pressure,
   * density, and composition are perfectly physical. Primitive recovery remains
   * the authoritative thermochemical consistency check after the update.
   */
  (void)starVelocity2;

  const su2double speciesTolerance = -1.0e-12 * starDensity;
  for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    if (starState[iSpecies] < speciesTolerance) return false;
  }

  return true;
}

CNumerics::ResidualType<> CUpwHLLC_Detonation::ComputeResidual(const CConfig*) {
  const su2double area = GeometryToolbox::Norm(nDim, Normal);
  if (area <= 0.0) {
    for (unsigned short iVar = 0; iVar < nVar; ++iVar) Flux[iVar] = 0.0;
    return ResidualType<>(Flux, nullptr, nullptr);
  }

  su2double unitNormal[MAXNDIM] = {0.0};
  for (unsigned short iDim = 0; iDim < nDim; ++iDim) unitNormal[iDim] = Normal[iDim] / area;

  GetInviscidProjFlux(U_i, V_i, unitNormal, FluxL);
  GetInviscidProjFlux(U_j, V_j, unitNormal, FluxR);

  const su2double rhoL = V_i[RHO_INDEX];
  const su2double rhoR = V_j[RHO_INDEX];
  const su2double pL = V_i[P_INDEX];
  const su2double pR = V_j[P_INDEX];
  const su2double aL = V_i[A_INDEX];
  const su2double aR = V_j[A_INDEX];
  const su2double unL = GeometryToolbox::DotProduct(nDim, &V_i[VEL_INDEX], unitNormal);
  const su2double unR = GeometryToolbox::DotProduct(nDim, &V_j[VEL_INDEX], unitNormal);

  MinPressure = std::min(MinPressure, std::min(pL, pR));
  MinDensity = std::min(MinDensity, std::min(rhoL, rhoR));

  const su2double sL = std::min(unL - aL, unR - aR);
  const su2double sR = std::max(unL + aL, unR + aR);
  const su2double denom = rhoL * (sL - unL) - rhoR * (sR - unR);

  if (rhoL <= 0.0 || rhoR <= 0.0 || pL <= 0.0 || pR <= 0.0 || aL <= 0.0 || aR <= 0.0 ||
      sL >= sR || fabs(denom) < 1.0e-14) {
    NegativeStateCount++;
    ComputeLLFFallback();
    return ResidualType<>(Flux, nullptr, nullptr);
  }

  const su2double sM = (pR - pL + rhoL * unL * (sL - unL) - rhoR * unR * (sR - unR)) / denom;

  su2double starPressureL = 0.0, starPressureR = 0.0;
  su2double starDensityL = 0.0, starDensityR = 0.0;

  if (0.0 <= sL) {
    for (unsigned short iVar = 0; iVar < nVar; ++iVar) Flux[iVar] = area * FluxL[iVar];
  } else if (sL <= 0.0 && 0.0 <= sM) {
    if (!ComputeStarState(U_i, V_i, unitNormal, sL, sM, unL, UStar, starPressureL, starDensityL)) {
      NegativeStateCount++;
      ComputeLLFFallback();
      return ResidualType<>(Flux, nullptr, nullptr);
    }
    for (unsigned short iVar = 0; iVar < nVar; ++iVar) {
      Flux[iVar] = area * (FluxL[iVar] + sL * (UStar[iVar] - U_i[iVar]));
    }
  } else if (sM <= 0.0 && 0.0 <= sR) {
    if (!ComputeStarState(U_j, V_j, unitNormal, sR, sM, unR, UStar, starPressureR, starDensityR)) {
      NegativeStateCount++;
      ComputeLLFFallback();
      return ResidualType<>(Flux, nullptr, nullptr);
    }
    for (unsigned short iVar = 0; iVar < nVar; ++iVar) {
      Flux[iVar] = area * (FluxR[iVar] + sR * (UStar[iVar] - U_j[iVar]));
    }
  } else {
    for (unsigned short iVar = 0; iVar < nVar; ++iVar) Flux[iVar] = area * FluxR[iVar];
  }

  return ResidualType<>(Flux, nullptr, nullptr);
}
