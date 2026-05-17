#pragma once

#include <array>
#include <memory>
#include <string>
#include <vector>

#include "CFVMFlowSolverBase.hpp"
#include "../fluid/CCVODEChemistryIntegrator.hpp"
#include "../fluid/CCanteraReactiveGas.hpp"
#include "../numerics/detonation/CDetonationNumerics.hpp"
#include "../variables/CDetonationEulerVariable.hpp"

class CDetonationEulerSolver final : public CFVMFlowSolverBase<CDetonationEulerVariable, ENUM_REGIME::COMPRESSIBLE> {
private:
  using BaseClass = CFVMFlowSolverBase<CDetonationEulerVariable, ENUM_REGIME::COMPRESSIBLE>;

  struct ChemEnergyAuditSample {
    su2double relError = 0.0;
    su2double clipCorrection = 0.0;
    su2double sumYCorrection = 0.0;
    su2double renormalizationCorrection = 0.0;
  };

  struct ChemEnergyWorstRecord {
    unsigned long timeIter = 0;
    su2double physicalTime = 0.0;
    unsigned long globalNode = 0;
    std::array<su2double, 3> coord{{0.0, 0.0, 0.0}};
    su2double density = 0.0;
    su2double temperatureBefore = 0.0;
    su2double temperatureAfter = 0.0;
    su2double pressureBefore = 0.0;
    su2double pressureAfter = 0.0;
    su2double internalEnergyBefore = 0.0;
    su2double internalEnergyAfter = 0.0;
    su2double chemistryEnergyAbsError = 0.0;
    su2double chemistryEnergyRelError = 0.0;
    su2double velocityMagnitude = 0.0;
    su2double sumYBefore = 0.0;
    su2double sumYAfter = 0.0;
    su2double minYBefore = 0.0;
    su2double minYAfter = 0.0;
    su2double maxClipCorrection = 0.0;
    su2double sumYCorrection = 0.0;
    su2double maxRenormalizationCorrection = 0.0;
    su2double renormalizationFactor = 1.0;
    su2double temperatureClipCorrection = 0.0;
    unsigned long retrySubsteps = 0;
    unsigned long cvodeInternalSteps = 0;
    int cvodeReturnFlag = 0;
    su2double qdot = 0.0;
    su2double maxAbsWdot = 0.0;
    std::array<int, 3> dominantSpecies{{-1, -1, -1}};
    std::array<su2double, 3> dominantWdot{{0.0, 0.0, 0.0}};
    std::vector<unsigned short> clippedSpecies;
  };

  unsigned short nSpecies = 0;
  CCanteraReactiveGas* ReactiveModel = nullptr;
  const CConfig* SolverConfig = nullptr;
  CNumerics** CachedNumerics = nullptr;
  std::unique_ptr<IChemistryIntegrator> ChemistryIntegrator;
  unsigned long ChemFailCells = 0;
  unsigned long ChemSubstepsMax = 0;
  su2double ChemEnergyAbsMax = 0.0;
  su2double ChemEnergyRelMax = 0.0;
  std::vector<su2double> ChemQdot;
  std::vector<unsigned long> ChemFailCount;
  std::vector<unsigned short> ChemFailFlag;
  std::vector<unsigned long> ChemSubstepsLast;
  std::vector<su2double> DebugSolutionSnapshot;
  std::vector<ChemEnergyAuditSample> ChemEnergyAuditSamples;
  std::vector<ChemEnergyWorstRecord> ChemEnergyWorstRecords;

  void SetNondimensionalization(CConfig* config);
  unsigned long SetPrimitive_Variables(CConfig* config);
  void AdvanceChemistry(CGeometry* geometry, CConfig* config);
  void AssembleFlowResidual(CGeometry* geometry, CSolver** solver_container, CConfig* config);
  void ApplyExplicitFlowStep(CGeometry* geometry, CSolver** solver_container, CConfig* config, su2double dtFactor);
  void RecomputeConservativeVector(const su2double* primitive, const su2double* massFractions, su2double* conservative) const;
  void BuildBoundaryStateFromPrimitive(const su2double* primitive, const su2double* massFractions, su2double* conservative) const;
  void ApplyVerificationLeftRightInitialization(CGeometry* geometry, CConfig* config);
  void SaveDebugSolutionSnapshot();
  void DebugResidualStage(const std::string& stage, CGeometry* geometry, CConfig* config,
                          const std::string& markerName = std::string()) const;
  void DebugPrimitiveStage(const std::string& stage, CGeometry* geometry, CConfig* config) const;
  void BuildSupersonicInletBoundaryState(CConfig* config, unsigned short marker,
                                         std::vector<su2double>& ghostPrimitive,
                                         std::vector<su2double>& ghostMassFractions,
                                         std::vector<su2double>& ghostConservative) const;
  inline bool DebugResidualBreakdownEnabled(const CConfig* config) const {
    return config->GetDetonation_Debug_Residual_Breakdown();
  }
  inline bool DebugBCResidualEnabled(const CConfig* config) const {
    return config->GetDetonation_Debug_BC_Residual();
  }
  inline bool ChemistryOnlyMode(const CConfig* config) const {
    return config->GetDetonation_Chemistry_Only();
  }
  inline bool ChemistryEnergyAuditEnabled(const CConfig* config) const {
    return config->GetDetonation_Chem_Energy_Audit();
  }
  inline bool UseStrangSplitting(const CConfig* config) const {
    return config->GetChemistry_Splitting() == "STRANG";
  }
  void WriteChemistryAuditFiles() const;

protected:
  void SetReferenceValues(const CConfig& config) final;
  void SetMax_Eigenvalue(CGeometry* geometry, const CConfig* config);

public:
  CDetonationEulerSolver() = delete;
  CDetonationEulerSolver(CGeometry* geometry, CConfig* config, unsigned short iMesh);
  ~CDetonationEulerSolver() override;

  void Preprocessing(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                     unsigned short iMesh, unsigned short iRKStep,
                     unsigned short RunTime_EqSystem, bool Output) final;

  void SetTime_Step(CGeometry* geometry, CSolver** solver_container,
                    CConfig* config, unsigned short iMesh, unsigned long Iteration) final;

  void Centered_Residual(CGeometry* geometry, CSolver** solver_container, CNumerics** numerics,
                         CConfig* config, unsigned short iMesh, unsigned short iRKStep) final;

  void Upwind_Residual(CGeometry* geometry, CSolver** solver_container, CNumerics** numerics,
                       CConfig* config, unsigned short iMesh) final;

  void Source_Residual(CGeometry* geometry, CSolver** solver_container, CNumerics** numerics,
                       CConfig* config, unsigned short iMesh) final;

  void ExplicitEuler_Iteration(CGeometry* geometry, CSolver** solver_container, CConfig* config) final;

  void BC_Supersonic_Inlet(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                           CNumerics* visc_numerics, CConfig* config, unsigned short marker) final;

  void BC_Supersonic_Outlet(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                            CNumerics* visc_numerics, CConfig* config, unsigned short marker) final;

  void BC_Far_Field(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                    CNumerics* visc_numerics, CConfig* config, unsigned short marker) final;

  void BC_Sym_Plane(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                    CNumerics* visc_numerics, CConfig* config, unsigned short marker) final;

  inline unsigned short GetNumSpecies() const { return nSpecies; }
  inline unsigned long GetChemFailCells() const { return ChemFailCells; }
  inline unsigned long GetChemSubstepsMax() const { return ChemSubstepsMax; }
  inline su2double GetChemEnergyAbsMax() const { return ChemEnergyAbsMax; }
  inline su2double GetChemEnergyRelMax() const { return ChemEnergyRelMax; }
  inline su2double GetChemQdot(unsigned long iPoint) const { return ChemQdot[iPoint]; }
  inline unsigned long GetChemFailCount(unsigned long iPoint) const { return ChemFailCount[iPoint]; }
  inline unsigned short GetChemFailFlag(unsigned long iPoint) const { return ChemFailFlag[iPoint]; }
  inline unsigned long GetChemSubstepsLast(unsigned long iPoint) const { return ChemSubstepsLast[iPoint]; }
  inline unsigned long GetHLLCFallbackCount() const { return CUpwHLLC_Detonation::GetFallbackCount(); }
  inline unsigned long GetHLLCNegativeStateCount() const { return CUpwHLLC_Detonation::GetNegativeStateCount(); }
  inline su2double GetHLLCMinPressure() const { return CUpwHLLC_Detonation::GetMinPressure(); }
  inline su2double GetHLLCMinDensity() const { return CUpwHLLC_Detonation::GetMinDensity(); }

  void PrintVerificationError(const CConfig*) const final {}
};
