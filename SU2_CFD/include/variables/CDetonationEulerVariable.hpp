#pragma once

#include <limits>
#include <vector>

#include "CFlowVariable.hpp"
#include "../fluid/CCanteraReactiveGas.hpp"

class CDetonationEulerVariable final : public CFlowVariable {
public:
  static constexpr size_t MAXNVAR = 128;

  template <class IndexType>
  struct CIndices {
    const IndexType nDim, nSpecies;
    CIndices(IndexType ndim, IndexType nspecies) : nDim(ndim), nSpecies(nspecies) {}
    inline IndexType NDim() const { return nDim; }
    inline IndexType NSpecies() const { return nSpecies; }
    inline IndexType SpeciesDensities() const { return 0; }
    inline IndexType Temperature() const { return nSpecies; }
    inline IndexType Velocity() const { return nSpecies + 1; }
    inline IndexType Pressure() const { return nSpecies + nDim + 1; }
    inline IndexType Density() const { return nSpecies + nDim + 2; }
    inline IndexType Enthalpy() const { return nSpecies + nDim + 3; }
    inline IndexType SoundSpeed() const { return nSpecies + nDim + 4; }
    inline IndexType LaminarViscosity() const { return nSpecies + nDim + 5; }
    inline IndexType EddyViscosity() const { return nSpecies + nDim + 6; }
    inline IndexType ThermalConductivity() const { return nSpecies + nDim + 7; }
    inline IndexType CpTotal() const { return nSpecies + nDim + 8; }
    inline IndexType Temperature_ve() const { return std::numeric_limits<IndexType>::max(); }
  };

private:
  const CIndices<unsigned long> indices;
  unsigned short nSpecies = 0;
  MatrixType Secondary;
  MatrixType MassFractions;
  VectorType MixtureGamma;

public:
  CDetonationEulerVariable(su2double pressure,
                           const su2double* massFractions,
                           const su2double* velocity,
                           su2double temperature,
                           unsigned long nPoint,
                           unsigned long nDim,
                           unsigned long nVar,
                           unsigned long nPrimVar,
                           unsigned long nPrimVarGrad,
                           const CConfig* config,
                           CCanteraReactiveGas* fluidModel);

  inline void SetdPdrho_e(unsigned long iPoint, su2double dPdrho_e) final { Secondary(iPoint, 0) = dPdrho_e; }
  inline void SetdPde_rho(unsigned long iPoint, su2double dPde_rho) final { Secondary(iPoint, 1) = dPde_rho; }

  inline bool SetDensity(unsigned long iPoint) final {
    su2double rho = 0.0;
    for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) rho += Solution(iPoint, iSpecies);
    Primitive(iPoint, indices.Density()) = rho;
    return rho <= 0.0;
  }

  inline bool SetPressure(unsigned long iPoint, su2double pressure) final {
    Primitive(iPoint, indices.Pressure()) = pressure;
    return pressure <= 0.0;
  }

  inline bool SetSoundSpeed(unsigned long iPoint, su2double soundspeed2) final {
    if (soundspeed2 < 0.0) return true;
    Primitive(iPoint, indices.SoundSpeed()) = sqrt(soundspeed2);
    return false;
  }

  inline bool SetTemperature(unsigned long iPoint, su2double temperature) final {
    Primitive(iPoint, indices.Temperature()) = temperature;
    return temperature <= 0.0;
  }

  inline void SetEnthalpy(unsigned long iPoint) final {
    Primitive(iPoint, indices.Enthalpy()) =
        (Solution(iPoint, nVar - 1) + Primitive(iPoint, indices.Pressure())) / Primitive(iPoint, indices.Density());
  }

  bool SetPrimVar(unsigned long iPoint, CFluidModel* fluidModel) final;

  void SetSecondaryVar(unsigned long iPoint, CFluidModel* fluidModel) override;

  inline su2double GetSecondary(unsigned long iPoint, unsigned long iVar) const final { return Secondary(iPoint, iVar); }
  inline su2double* GetSecondary(unsigned long iPoint) final { return Secondary[iPoint]; }
  inline void SetSecondary(unsigned long iPoint, unsigned long iVar, su2double value) final { Secondary(iPoint, iVar) = value; }
  inline void SetSecondary(unsigned long iPoint, const su2double* secondary) final {
    for (unsigned long iVar = 0; iVar < nSecondaryVar; ++iVar) Secondary(iPoint, iVar) = secondary[iVar];
  }

  inline void SetVelocity(unsigned long iPoint) final {
    Velocity2(iPoint) = 0.0;
    const su2double rho = Primitive(iPoint, indices.Density());
    for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
      Primitive(iPoint, indices.Velocity() + iDim) = Solution(iPoint, nSpecies + iDim) / rho;
      Velocity2(iPoint) += Primitive(iPoint, indices.Velocity() + iDim) * Primitive(iPoint, indices.Velocity() + iDim);
    }
  }

  inline su2double GetPressure(unsigned long iPoint) const final { return Primitive(iPoint, indices.Pressure()); }
  inline su2double GetSoundSpeed(unsigned long iPoint) const final { return Primitive(iPoint, indices.SoundSpeed()); }
  inline su2double GetEnthalpy(unsigned long iPoint) const final { return Primitive(iPoint, indices.Enthalpy()); }
  inline su2double GetDensity(unsigned long iPoint) const final { return Primitive(iPoint, indices.Density()); }
  inline su2double GetEnergy(unsigned long iPoint) const final { return Solution(iPoint, nVar - 1) / Primitive(iPoint, indices.Density()); }
  inline su2double GetTemperature(unsigned long iPoint) const final { return Primitive(iPoint, indices.Temperature()); }
  inline su2double GetVelocity(unsigned long iPoint, unsigned long iDim) const final {
    return Primitive(iPoint, indices.Velocity() + iDim);
  }

  inline su2double GetMassFraction(unsigned long iPoint, unsigned short iSpecies) const {
    return MassFractions(iPoint, iSpecies);
  }

  inline const su2double* GetMassFractions(unsigned long iPoint) const { return MassFractions[iPoint]; }

  inline su2double GetGamma(unsigned long iPoint) const { return MixtureGamma(iPoint); }

  inline unsigned short GetNumSpecies() const { return nSpecies; }
};
