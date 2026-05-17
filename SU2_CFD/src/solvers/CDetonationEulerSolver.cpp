#include "../../include/solvers/CDetonationEulerSolver.hpp"

#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>

#include "../../include/numerics/detonation/CDetonationNumerics.hpp"
#include "../../include/solvers/CFVMFlowSolverBase.inl"
#include "../../../Common/include/option_structure.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"

template class CFVMFlowSolverBase<CDetonationEulerVariable, ENUM_REGIME::COMPRESSIBLE>;

namespace {

std::string BoundaryKindName(unsigned short kind) {
  switch (kind) {
    case EULER_WALL: return "EULER_WALL";
    case FAR_FIELD: return "FAR_FIELD";
    case SYMMETRY_PLANE: return "SYMMETRY_PLANE";
    case SUPERSONIC_INLET: return "SUPERSONIC_INLET";
    case SUPERSONIC_OUTLET: return "SUPERSONIC_OUTLET";
    default: return "OTHER";
  }
}

std::string EscapeCSV(const std::string& value) {
  std::string escaped;
  escaped.reserve(value.size());
  for (char c : value) {
    if (c == '"') escaped.push_back('"');
    escaped.push_back(c);
  }
  return escaped;
}

std::string SerializeArray(const su2double* data, unsigned short count) {
  std::ostringstream out;
  out << std::scientific << std::setprecision(17);
  out << '[';
  for (unsigned short i = 0; i < count; ++i) {
    if (i != 0) out << ';';
    out << data[i];
  }
  out << ']';
  return out.str();
}

std::string SerializeIndexVector(const std::vector<unsigned short>& values) {
  std::ostringstream out;
  out << '[';
  for (unsigned long i = 0; i < values.size(); ++i) {
    if (i != 0) out << ';';
    out << values[i];
  }
  out << ']';
  return out.str();
}

std::string InsertSuffixBeforeExtension(const std::string& path, const std::string& suffix) {
  const auto pos = path.find_last_of('.');
  if (pos == std::string::npos || pos == 0 || path.find('/') > pos) return path + suffix;
  return path.substr(0, pos) + suffix + path.substr(pos);
}

void ComputeProjectedFlux(unsigned short nSpecies, unsigned short nDim, unsigned short nVar,
                          const su2double* conserved, const su2double* primitive,
                          const su2double* normal, su2double* projectedFlux) {
  const unsigned short rhoIndex = nSpecies + nDim + 2;
  const unsigned short velIndex = nSpecies + 1;
  const unsigned short pressureIndex = nSpecies + nDim + 1;
  const unsigned short enthalpyIndex = nSpecies + nDim + 3;

  for (unsigned short iVar = 0; iVar < nVar; ++iVar) projectedFlux[iVar] = 0.0;

  su2double projVelocity = 0.0;
  for (unsigned short iDim = 0; iDim < nDim; ++iDim) projVelocity += primitive[velIndex + iDim] * normal[iDim];

  for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    projectedFlux[iSpecies] = conserved[iSpecies] * projVelocity;
  }

  for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
    projectedFlux[nSpecies + iDim] =
        conserved[nSpecies + iDim] * projVelocity + primitive[pressureIndex] * normal[iDim];
  }

  projectedFlux[nVar - 1] = primitive[rhoIndex] * projVelocity * primitive[enthalpyIndex];
}

struct BoundaryDebugSummary {
  bool enabled = false;
  unsigned short nSpecies = 0;
  unsigned short nDim = 0;
  unsigned short nVar = 0;
  unsigned short nPrimVar = 0;
  unsigned long localFaceCount = 0;
  su2double maxSpeciesAbs = -1.0;
  su2double maxEnergyAbs = -1.0;
  su2double maxTraceAbs = -1.0;
  su2double traceArea = 0.0;
  su2double traceUNInterior = 0.0;
  su2double traceUNExterior = 0.0;
  std::vector<su2double> maxMomentumAbs;
  std::vector<su2double> sumResidual;
  std::vector<unsigned long> maxMomentumPoint;
  std::vector<unsigned long> maxMomentumFace;
  unsigned long maxSpeciesPoint = 0;
  unsigned long maxEnergyPoint = 0;
  unsigned long maxTracePoint = 0;
  unsigned long maxSpeciesFace = 0;
  unsigned long maxEnergyFace = 0;
  unsigned long maxTraceFace = 0;
  std::vector<su2double> maxSpeciesCoord;
  std::vector<su2double> maxEnergyCoord;
  std::vector<su2double> maxTraceCoord;
  std::vector<su2double> traceNormal;
  std::vector<su2double> tracePrimitiveInterior;
  std::vector<su2double> tracePrimitiveExterior;
  std::vector<su2double> traceFluxInterior;
  std::vector<su2double> traceFluxExterior;
  std::vector<su2double> traceResidual;

  BoundaryDebugSummary(bool enabledIn, unsigned short nSpeciesIn, unsigned short nDimIn,
                       unsigned short nVarIn, unsigned short nPrimVarIn)
      : enabled(enabledIn),
        nSpecies(nSpeciesIn),
        nDim(nDimIn),
        nVar(nVarIn),
        nPrimVar(nPrimVarIn),
        maxMomentumAbs(nDimIn, -1.0),
        sumResidual(nVarIn, 0.0),
        maxMomentumPoint(nDimIn, 0),
        maxMomentumFace(nDimIn, 0),
        maxSpeciesCoord(nDimIn, 0.0),
        maxEnergyCoord(nDimIn, 0.0),
        maxTraceCoord(nDimIn, 0.0),
        traceNormal(nDimIn, 0.0),
        tracePrimitiveInterior(nPrimVarIn, 0.0),
        tracePrimitiveExterior(nPrimVarIn, 0.0),
        traceFluxInterior(nVarIn, 0.0),
        traceFluxExterior(nVarIn, 0.0),
        traceResidual(nVarIn, 0.0) {}

  void Accumulate(unsigned long localFace, unsigned long localPoint, unsigned long globalPoint, const su2double* coord,
                  const su2double* normal, const su2double* primitiveInterior,
                  const su2double* primitiveExterior, const su2double* conservativeInterior,
                  const su2double* conservativeExterior, const CNumerics::ResidualType<>& residual) {
    if (!enabled) return;

    localFaceCount++;

    for (unsigned short iVar = 0; iVar < nVar; ++iVar) sumResidual[iVar] += residual[iVar];

    for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      const su2double value = fabs(residual[iSpecies]);
      if (value > maxSpeciesAbs) {
        maxSpeciesAbs = value;
        maxSpeciesPoint = globalPoint;
        maxSpeciesFace = localFace;
        for (unsigned short iDim = 0; iDim < nDim; ++iDim) maxSpeciesCoord[iDim] = coord[iDim];
      }
    }

    for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
      const su2double value = fabs(residual[nSpecies + iDim]);
      if (value > maxMomentumAbs[iDim]) {
        maxMomentumAbs[iDim] = value;
        maxMomentumPoint[iDim] = globalPoint;
        maxMomentumFace[iDim] = localFace;
      }
    }

    const su2double energyValue = fabs(residual[nVar - 1]);
    if (energyValue > maxEnergyAbs) {
      maxEnergyAbs = energyValue;
      maxEnergyPoint = globalPoint;
      maxEnergyFace = localFace;
      for (unsigned short iDim = 0; iDim < nDim; ++iDim) maxEnergyCoord[iDim] = coord[iDim];
    }

    su2double traceValue = 0.0;
    for (unsigned short iVar = 0; iVar < nVar; ++iVar) traceValue = max(traceValue, fabs(residual[iVar]));

    if (traceValue > maxTraceAbs) {
      maxTraceAbs = traceValue;
      maxTracePoint = globalPoint;
      maxTraceFace = localFace;
      for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
        maxTraceCoord[iDim] = coord[iDim];
        traceNormal[iDim] = normal[iDim];
      }
      traceArea = GeometryToolbox::Norm(nDim, normal);
      traceUNInterior = 0.0;
      traceUNExterior = 0.0;
      const unsigned short velIndex = nSpecies + 1;
      for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
        traceUNInterior += primitiveInterior[velIndex + iDim] * normal[iDim];
        traceUNExterior += primitiveExterior[velIndex + iDim] * normal[iDim];
      }
      for (unsigned short iPrim = 0; iPrim < nPrimVar; ++iPrim) {
        tracePrimitiveInterior[iPrim] = primitiveInterior[iPrim];
        tracePrimitiveExterior[iPrim] = primitiveExterior[iPrim];
      }
      for (unsigned short iVar = 0; iVar < nVar; ++iVar) traceResidual[iVar] = residual[iVar];
      ComputeProjectedFlux(nSpecies, nDim, nVar, conservativeInterior, primitiveInterior, normal, traceFluxInterior.data());
      ComputeProjectedFlux(nSpecies, nDim, nVar, conservativeExterior, primitiveExterior, normal, traceFluxExterior.data());
    }
  }
};

void InitializeBCResidualCSV(const CConfig* config) {
  static bool initialized = false;
  if (!config->GetDetonation_Debug_BC_Residual()) return;
  if (initialized) return;

  const int rank = SU2_MPI::GetRank();
  if (rank == MASTER_NODE) {
    std::ofstream out("bc_residual_breakdown.csv", std::ios::trunc);
    out << "time_iter,physical_time,rank,marker,marker_kind,local_boundary_faces,"
        << "max_abs_species_res,max_abs_species_node_global,max_abs_species_face_local,max_abs_species_coord,"
        << "max_abs_xmom_res,max_abs_xmom_node_global,max_abs_xmom_face_local,"
        << "max_abs_ymom_res,max_abs_ymom_node_global,max_abs_ymom_face_local,"
        << "max_abs_zmom_res,max_abs_zmom_node_global,max_abs_zmom_face_local,"
        << "max_abs_energy_res,max_abs_energy_node_global,max_abs_energy_face_local,max_abs_energy_coord,"
        << "sum_species_res,sum_xmom_res,sum_ymom_res,sum_zmom_res,sum_energy_res,"
        << "trace_node_global,trace_face_local,trace_coord,trace_normal,trace_area,trace_u_n_interior,trace_u_n_exterior,"
        << "trace_primitive_interior,trace_primitive_exterior,trace_physical_flux_interior,trace_physical_flux_exterior,"
        << "trace_residual_added\n";
  }
  SU2_MPI::Barrier(SU2_MPI::GetComm());
  initialized = true;
}

void AppendBCResidualCSV(const CGeometry* geometry, const CConfig* config, unsigned short marker,
                         const BoundaryDebugSummary& summary) {
  if (!summary.enabled) return;

  InitializeBCResidualCSV(config);

  const int rank = SU2_MPI::GetRank();
  int size = 1;
#ifdef HAVE_MPI
  SU2_MPI::Comm_size(SU2_MPI::GetComm(), &size);
#endif

  const std::string markerName = config->GetMarker_All_TagBound(marker);
  const std::string markerKind = BoundaryKindName(config->GetMarker_All_KindBC(marker));

  const auto coordString = [&](const std::vector<su2double>& coord) {
    return SerializeArray(coord.data(), static_cast<unsigned short>(coord.size()));
  };

  for (int rankPrint = 0; rankPrint < size; ++rankPrint) {
    if (rank == rankPrint) {
      std::ofstream out("bc_residual_breakdown.csv", std::ios::app);
      out << std::scientific << std::setprecision(17);
      out << config->GetTimeIter() << ','
          << config->GetPhysicalTime() << ','
          << rank << ','
          << '"' << EscapeCSV(markerName) << '"' << ','
          << '"' << EscapeCSV(markerKind) << '"' << ','
          << summary.localFaceCount << ','
          << summary.maxSpeciesAbs << ','
          << summary.maxSpeciesPoint << ','
          << summary.maxSpeciesFace << ','
          << '"' << EscapeCSV(coordString(summary.maxSpeciesCoord)) << '"' << ','
          << ((summary.nDim > 0) ? summary.maxMomentumAbs[0] : 0.0) << ','
          << ((summary.nDim > 0) ? summary.maxMomentumPoint[0] : 0ul) << ','
          << ((summary.nDim > 0) ? summary.maxMomentumFace[0] : 0ul) << ','
          << ((summary.nDim > 1) ? summary.maxMomentumAbs[1] : 0.0) << ','
          << ((summary.nDim > 1) ? summary.maxMomentumPoint[1] : 0ul) << ','
          << ((summary.nDim > 1) ? summary.maxMomentumFace[1] : 0ul) << ','
          << ((summary.nDim > 2) ? summary.maxMomentumAbs[2] : 0.0) << ','
          << ((summary.nDim > 2) ? summary.maxMomentumPoint[2] : 0ul) << ','
          << ((summary.nDim > 2) ? summary.maxMomentumFace[2] : 0ul) << ','
          << summary.maxEnergyAbs << ','
          << summary.maxEnergyPoint << ','
          << summary.maxEnergyFace << ','
          << '"' << EscapeCSV(coordString(summary.maxEnergyCoord)) << '"' << ',';

      su2double sumSpecies = 0.0;
      for (unsigned short iSpecies = 0; iSpecies < summary.nSpecies; ++iSpecies) sumSpecies += summary.sumResidual[iSpecies];
      out << sumSpecies << ','
          << ((summary.nDim > 0) ? summary.sumResidual[summary.nSpecies] : 0.0) << ','
          << ((summary.nDim > 1) ? summary.sumResidual[summary.nSpecies + 1] : 0.0) << ','
          << ((summary.nDim > 2) ? summary.sumResidual[summary.nSpecies + 2] : 0.0) << ','
          << summary.sumResidual[summary.nVar - 1] << ','
          << summary.maxTracePoint << ','
          << summary.maxTraceFace << ','
          << '"' << EscapeCSV(coordString(summary.maxTraceCoord)) << '"' << ','
          << '"' << EscapeCSV(coordString(summary.traceNormal)) << '"' << ','
          << summary.traceArea << ','
          << summary.traceUNInterior << ','
          << summary.traceUNExterior << ','
          << '"' << EscapeCSV(SerializeArray(summary.tracePrimitiveInterior.data(), summary.nPrimVar)) << '"' << ','
          << '"' << EscapeCSV(SerializeArray(summary.tracePrimitiveExterior.data(), summary.nPrimVar)) << '"' << ','
          << '"' << EscapeCSV(SerializeArray(summary.traceFluxInterior.data(), summary.nVar)) << '"' << ','
          << '"' << EscapeCSV(SerializeArray(summary.traceFluxExterior.data(), summary.nVar)) << '"' << ','
          << '"' << EscapeCSV(SerializeArray(summary.traceResidual.data(), summary.nVar)) << '"'
          << '\n';
    }
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }
}

}  // namespace

CDetonationEulerSolver::CDetonationEulerSolver(CGeometry* geometry, CConfig* config, unsigned short iMesh)
    : CFVMFlowSolverBase<CDetonationEulerVariable, ENUM_REGIME::COMPRESSIBLE>(*geometry, *config) {
  SolverConfig = config;
  MGLevel = iMesh;
  dynamic_grid = config->GetDynamic_Grid();
  nSpecies = config->GetnSpecies();
  nMarker = config->GetnMarker_All();
  nDim = geometry->GetnDim();
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  nVar = nSpecies + nDim + 1;
  nPrimVar = nSpecies + nDim + 5;
  nPrimVarGrad = nPrimVar;
  nVarGrad = nPrimVarGrad;

  nVertex.resize(nMarker);
  for (auto iMarker = 0ul; iMarker < nMarker; ++iMarker) nVertex[iMarker] = geometry->nVertex[iMarker];

  ReactiveModel = new CCanteraReactiveGas(config, nDim);
  SetNondimensionalization(config);
  ChemistryIntegrator = std::unique_ptr<IChemistryIntegrator>(new CCVODEChemistryIntegrator(*ReactiveModel, config));
  ChemQdot.assign(nPoint, 0.0);
  ChemFailCount.assign(nPoint, 0ul);
  ChemFailFlag.assign(nPoint, 0u);
  ChemSubstepsLast.assign(nPoint, 0ul);
  DebugSolutionSnapshot.assign(nPointDomain * nVar, 0.0);

  AllocateTerribleLegacyTemporaryVariables();
  Allocate(*config);

  const su2double* massFractions = config->GetGas_Composition();
  nodes = new CDetonationEulerVariable(config->GetPressure_FreeStreamND(),
                                       massFractions,
                                       config->GetVelocity_FreeStreamND(),
                                       config->GetTemperature_FreeStreamND(),
                                       nPoint, nDim, nVar, nPrimVar, nPrimVarGrad, config, ReactiveModel);
  SetBaseClassPointerToNodes();
  ApplyVerificationLeftRightInitialization(geometry, config);
  CommunicateInitialState(geometry, config);
  SolverName = "DETONATION.FLOW";

  if (ChemistryOnlyMode(config) && SU2_MPI::GetRank() == MASTER_NODE) {
    cout << "WARNING: DETONATION_CHEMISTRY_ONLY is for 0D verification only." << endl;
  }
}

CDetonationEulerSolver::~CDetonationEulerSolver() {
  WriteChemistryAuditFiles();
  delete ReactiveModel;
}

void CDetonationEulerSolver::ApplyVerificationLeftRightInitialization(CGeometry* geometry, CConfig* config) {
  if (!config->GetDetonation_LeftRight_Init()) return;

  std::vector<su2double> leftPrimitive(nPrimVar, 0.0), rightPrimitive(nPrimVar, 0.0);
  std::vector<su2double> leftConservative(nVar, 0.0), rightConservative(nVar, 0.0);
  std::vector<su2double> massFractions(nSpecies, 0.0);
  for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    massFractions[iSpecies] = config->GetGas_Composition()[iSpecies];
  }

  auto fillPrimitive = [&](su2double pressure, su2double temperature,
                           su2double velocityX, su2double velocityY, su2double velocityZ,
                           std::vector<su2double>& primitive, std::vector<su2double>& conservative) {
    ReactiveModel->SetMassFractions(massFractions.data(), nSpecies);
    ReactiveModel->SetTDState_PT(pressure, temperature);

    const su2double density = ReactiveModel->GetDensity();
    const su2double velocity[3] = {velocityX, velocityY, velocityZ};
    su2double velocity2 = 0.0;
    for (unsigned short iDim = 0; iDim < nDim; ++iDim) velocity2 += velocity[iDim] * velocity[iDim];

    for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      primitive[iSpecies] = density * massFractions[iSpecies];
    }
    primitive[nSpecies] = temperature;
    for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
      primitive[nSpecies + 1 + iDim] = velocity[iDim];
    }
    primitive[nSpecies + nDim + 1] = pressure;
    primitive[nSpecies + nDim + 2] = density;
    primitive[nSpecies + nDim + 3] =
        ReactiveModel->GetStaticEnergy() + 0.5 * velocity2 + pressure / density;
    primitive[nSpecies + nDim + 4] = ReactiveModel->GetSoundSpeed();

    RecomputeConservativeVector(primitive.data(), massFractions.data(), conservative.data());
  };

  fillPrimitive(config->GetDetonation_Init_Left_Pressure(),
                config->GetDetonation_Init_Left_Temperature(),
                config->GetDetonation_Init_Left_U(),
                config->GetDetonation_Init_Left_V(),
                config->GetDetonation_Init_Left_W(),
                leftPrimitive, leftConservative);
  fillPrimitive(config->GetDetonation_Init_Right_Pressure(),
                config->GetDetonation_Init_Right_Temperature(),
                config->GetDetonation_Init_Right_U(),
                config->GetDetonation_Init_Right_V(),
                config->GetDetonation_Init_Right_W(),
                rightPrimitive, rightConservative);

  const su2double discontinuityX = config->GetDetonation_Init_Discontinuity_X();
  for (unsigned long iPoint = 0; iPoint < nPoint; ++iPoint) {
    if (!geometry->nodes->GetDomain(iPoint)) continue;
    const bool useLeft = geometry->nodes->GetCoord(iPoint, 0) < discontinuityX;
    const auto& primitive = useLeft ? leftPrimitive : rightPrimitive;
    const auto& conservative = useLeft ? leftConservative : rightConservative;
    nodes->SetSolution(iPoint, conservative.data());
    nodes->SetSolution_Old(iPoint, conservative.data());
    for (unsigned long iVar = 0; iVar < nPrimVar; ++iVar) {
      nodes->SetPrimitive(iPoint, iVar, primitive[iVar]);
    }
  }

  if (SU2_MPI::GetRank() == MASTER_NODE) {
    cout << "WARNING: DETONATION_LEFT_RIGHT_INIT is a verification-only initialization path." << endl;
  }
}

void CDetonationEulerSolver::SetReferenceValues(const CConfig& config) {
  DynamicPressureRef = 0.5 * config.GetDensity_FreeStreamND() *
                       config.GetModVel_FreeStreamND() * config.GetModVel_FreeStreamND();
  AeroCoeffForceRef = DynamicPressureRef * config.GetRefArea();
}

void CDetonationEulerSolver::SetNondimensionalization(CConfig* config) {
  const su2double pressure = config->GetPressure_FreeStream();
  const su2double temperature = config->GetTemperature_FreeStream();
  ReactiveModel->SetMassFractions(config->GetGas_Composition(), nSpecies);
  ReactiveModel->SetTDState_PT(pressure, temperature);

  config->SetDensity_FreeStream(ReactiveModel->GetDensity());
  const su2double soundSpeed = ReactiveModel->GetSoundSpeed();
  const su2double alpha = config->GetAoA() * PI_NUMBER / 180.0;
  const su2double beta = config->GetAoS() * PI_NUMBER / 180.0;
  const su2double mach = config->GetMach();

  if (nDim == 2) {
    config->GetVelocity_FreeStream()[0] = cos(alpha) * mach * soundSpeed;
    config->GetVelocity_FreeStream()[1] = sin(alpha) * mach * soundSpeed;
  } else {
    config->GetVelocity_FreeStream()[0] = cos(alpha) * cos(beta) * mach * soundSpeed;
    config->GetVelocity_FreeStream()[1] = sin(beta) * mach * soundSpeed;
    config->GetVelocity_FreeStream()[2] = sin(alpha) * cos(beta) * mach * soundSpeed;
  }

  su2double modVel = 0.0;
  for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
    modVel += config->GetVelocity_FreeStream()[iDim] * config->GetVelocity_FreeStream()[iDim];
  }
  modVel = sqrt(modVel);
  config->SetModVel_FreeStream(modVel);

  config->SetPressure_Ref(1.0);
  config->SetDensity_Ref(1.0);
  config->SetTemperature_Ref(1.0);
  config->SetVelocity_Ref(1.0);
  config->SetTime_Ref(1.0);
  config->SetGas_Constant_Ref(1.0);
  config->SetGas_ConstantND(config->GetGas_Constant());
  config->SetPressure_FreeStreamND(config->GetPressure_FreeStream());
  config->SetDensity_FreeStreamND(config->GetDensity_FreeStream());
  config->SetTemperature_FreeStreamND(config->GetTemperature_FreeStream());
  for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
    config->SetVelocity_FreeStreamND(config->GetVelocity_FreeStream()[iDim], iDim);
  }
  config->SetModVel_FreeStreamND(modVel);
}

void CDetonationEulerSolver::SaveDebugSolutionSnapshot() {
  for (unsigned long iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    for (unsigned short iVar = 0; iVar < nVar; ++iVar) {
      DebugSolutionSnapshot[iPoint * nVar + iVar] = nodes->GetSolution(iPoint, iVar);
    }
  }
}

void CDetonationEulerSolver::DebugResidualStage(const std::string& stage, CGeometry* geometry, CConfig* config,
                                                const std::string& markerName) const {
  if (!DebugResidualBreakdownEnabled(config)) return;

  const bool useActualDelta = (stage.find("after explicit update") != std::string::npos) ||
                              (stage.find("after chemistry") != std::string::npos) ||
                              (stage.find("after SetPrimitive") != std::string::npos);

  su2double speciesResMax = -1.0, energyResMax = -1.0, deltaUMax = -1.0;
  std::vector<su2double> momentumResMax(nDim, -1.0);
  unsigned long speciesPoint = 0, energyPoint = 0, deltaPoint = 0;
  std::vector<unsigned long> momentumPoint(nDim, 0);

  for (unsigned long iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    const su2double* residual = LinSysRes.GetBlock(iPoint);
    const su2double vol = geometry->nodes->GetVolume(iPoint) + geometry->nodes->GetPeriodicVolume(iPoint);
    const su2double deltaScale = (vol > 0.0) ? nodes->GetDelta_Time(iPoint) / vol : 0.0;

    for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      const su2double value = fabs(residual[iSpecies]);
      if (value > speciesResMax) {
        speciesResMax = value;
        speciesPoint = iPoint;
      }
    }

    for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
      const su2double value = fabs(residual[nSpecies + iDim]);
      if (value > momentumResMax[iDim]) {
        momentumResMax[iDim] = value;
        momentumPoint[iDim] = iPoint;
      }
    }

    const su2double energyValue = fabs(residual[nVar - 1]);
    if (energyValue > energyResMax) {
      energyResMax = energyValue;
      energyPoint = iPoint;
    }

    for (unsigned short iVar = 0; iVar < nVar; ++iVar) {
      su2double deltaU = 0.0;
      if (useActualDelta) {
        deltaU = fabs(nodes->GetSolution(iPoint, iVar) - DebugSolutionSnapshot[iPoint * nVar + iVar]);
      } else {
        deltaU = fabs(residual[iVar] * deltaScale);
      }
      if (deltaU > deltaUMax) {
        deltaUMax = deltaU;
        deltaPoint = iPoint;
      }
    }
  }

  const int rank = SU2_MPI::GetRank();
  int size = 1;
#ifdef HAVE_MPI
  SU2_MPI::Comm_size(SU2_MPI::GetComm(), &size);
#endif

  for (int rankPrint = 0; rankPrint < size; ++rankPrint) {
    if (rank == rankPrint) {
      std::ostringstream out;
      out << std::scientific << std::setprecision(12);
      out << "[DETONATION_DEBUG] rank=" << rank << " stage=\"" << stage << "\"";
      if (!markerName.empty()) out << " marker=\"" << markerName << "\"";
      out << '\n';

      auto appendPoint = [&](const std::string& name, su2double value, unsigned long point) {
        const auto coord = geometry->nodes->GetCoord(point);
        out << "  " << name
            << " value=" << value
            << " global_node=" << geometry->nodes->GetGlobalIndex(point)
            << " local_node=" << point
            << " coord=(" << coord[0];
        if (nDim > 1) out << ", " << coord[1];
        if (nDim > 2) out << ", " << coord[2];
        out << ")\n";
      };

      appendPoint("max|Res[rho_s]|", speciesResMax, speciesPoint);
      for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
        appendPoint("max|Res[rho*u_" + std::to_string(iDim) + "]|", momentumResMax[iDim], momentumPoint[iDim]);
      }
      appendPoint("max|Res[rhoE]|", energyResMax, energyPoint);
      appendPoint("max|deltaU|", deltaUMax, deltaPoint);
      cout << out.str() << std::flush;
    }
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }
}

void CDetonationEulerSolver::DebugPrimitiveStage(const std::string& stage, CGeometry* geometry, CConfig* config) const {
  if (!DebugResidualBreakdownEnabled(config)) return;

  su2double minRho = std::numeric_limits<su2double>::max();
  su2double maxRho = -std::numeric_limits<su2double>::max();
  su2double minSumRhoS = std::numeric_limits<su2double>::max();
  su2double maxSumRhoS = -std::numeric_limits<su2double>::max();
  su2double minRhoE = std::numeric_limits<su2double>::max();
  su2double maxRhoE = -std::numeric_limits<su2double>::max();
  su2double minE = std::numeric_limits<su2double>::max();
  su2double maxE = -std::numeric_limits<su2double>::max();
  su2double minT = std::numeric_limits<su2double>::max();
  su2double maxT = -std::numeric_limits<su2double>::max();
  su2double minP = std::numeric_limits<su2double>::max();
  su2double maxP = -std::numeric_limits<su2double>::max();
  su2double minA = std::numeric_limits<su2double>::max();
  su2double maxA = -std::numeric_limits<su2double>::max();
  su2double minH = std::numeric_limits<su2double>::max();
  su2double maxH = -std::numeric_limits<su2double>::max();
  su2double minSumYError = std::numeric_limits<su2double>::max();
  su2double maxSumYError = -std::numeric_limits<su2double>::max();
  std::vector<su2double> minVel(nDim, std::numeric_limits<su2double>::max());
  std::vector<su2double> maxVel(nDim, -std::numeric_limits<su2double>::max());

  for (unsigned long iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    su2double sumRhoS = 0.0, sumY = 0.0, vel2 = 0.0;
    for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      sumRhoS += nodes->GetSolution(iPoint, iSpecies);
      sumY += nodes->GetMassFraction(iPoint, iSpecies);
    }
    for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
      const su2double vel = nodes->GetVelocity(iPoint, iDim);
      minVel[iDim] = min(minVel[iDim], vel);
      maxVel[iDim] = max(maxVel[iDim], vel);
      vel2 += vel * vel;
    }

    const su2double rho = nodes->GetDensity(iPoint);
    const su2double rhoE = nodes->GetSolution(iPoint, nVar - 1);
    const su2double e = rhoE / rho - 0.5 * vel2;

    minRho = min(minRho, rho);
    maxRho = max(maxRho, rho);
    minSumRhoS = min(minSumRhoS, sumRhoS);
    maxSumRhoS = max(maxSumRhoS, sumRhoS);
    minRhoE = min(minRhoE, rhoE);
    maxRhoE = max(maxRhoE, rhoE);
    minE = min(minE, e);
    maxE = max(maxE, e);
    minT = min(minT, nodes->GetTemperature(iPoint));
    maxT = max(maxT, nodes->GetTemperature(iPoint));
    minP = min(minP, nodes->GetPressure(iPoint));
    maxP = max(maxP, nodes->GetPressure(iPoint));
    minA = min(minA, nodes->GetSoundSpeed(iPoint));
    maxA = max(maxA, nodes->GetSoundSpeed(iPoint));
    minH = min(minH, nodes->GetEnthalpy(iPoint));
    maxH = max(maxH, nodes->GetEnthalpy(iPoint));
    minSumYError = min(minSumYError, fabs(sumY - 1.0));
    maxSumYError = max(maxSumYError, fabs(sumY - 1.0));
  }

  const int rank = SU2_MPI::GetRank();
  int size = 1;
#ifdef HAVE_MPI
  SU2_MPI::Comm_size(SU2_MPI::GetComm(), &size);
#endif

  for (int rankPrint = 0; rankPrint < size; ++rankPrint) {
    if (rank == rankPrint) {
      std::ostringstream out;
      out << std::scientific << std::setprecision(12);
      out << "[DETONATION_DEBUG] rank=" << rank << " primitive-stage=\"" << stage << "\"\n";
      out << "  span(sum rho_s)=" << (maxSumRhoS - minSumRhoS) << '\n';
      out << "  span(rho)=" << (maxRho - minRho) << '\n';
      out << "  span(rhoE)=" << (maxRhoE - minRhoE) << '\n';
      out << "  span(e)=" << (maxE - minE) << '\n';
      out << "  span(T)=" << (maxT - minT) << '\n';
      out << "  span(p)=" << (maxP - minP) << '\n';
      out << "  span(a)=" << (maxA - minA) << '\n';
      out << "  span(h)=" << (maxH - minH) << '\n';
      out << "  max|sum(Y)-1|=" << maxSumYError << '\n';
      for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
        out << "  span(u_" << iDim << ")=" << (maxVel[iDim] - minVel[iDim]) << '\n';
      }
      cout << out.str() << std::flush;
    }
    SU2_MPI::Barrier(SU2_MPI::GetComm());
  }
}

void CDetonationEulerSolver::WriteChemistryAuditFiles() const {
  if (SolverConfig == nullptr || !ChemistryEnergyAuditEnabled(SolverConfig)) return;

  const std::string basePath = SolverConfig->GetDetonation_Chem_Energy_Audit_File();
  if (basePath.empty() || basePath == "NONE") return;

  const int rank = SU2_MPI::GetRank();
  const std::string rankSuffix = ".rank" + std::to_string(rank);
  const std::string worstPath = InsertSuffixBeforeExtension(basePath, rankSuffix);
  const std::string samplesPath = InsertSuffixBeforeExtension(basePath, ".all_cells" + rankSuffix);

  {
    std::ofstream out(worstPath);
    if (out) {
      out << "time_iter,physical_time,global_node,x,y,z,rho,T_before,T_after,p_before,p_after,"
             "e_before,e_after,e_abs_error,e_rel_error,velocity_mag,sumY_before,sumY_after,"
             "minY_before,minY_after,max_clip_correction,sumY_correction,max_renorm_correction,"
             "renorm_factor,temp_clip_correction,retry_substeps,cvode_internal_steps,cvode_return_flag,"
             "qdot,max_abs_wdot,dominant_species_0,dominant_wdot_0,dominant_species_1,dominant_wdot_1,"
             "dominant_species_2,dominant_wdot_2,clipped_species\n";
      out << std::scientific << std::setprecision(17);
      for (const auto& record : ChemEnergyWorstRecords) {
        out << record.timeIter << ','
            << record.physicalTime << ','
            << record.globalNode << ','
            << record.coord[0] << ','
            << record.coord[1] << ','
            << record.coord[2] << ','
            << record.density << ','
            << record.temperatureBefore << ','
            << record.temperatureAfter << ','
            << record.pressureBefore << ','
            << record.pressureAfter << ','
            << record.internalEnergyBefore << ','
            << record.internalEnergyAfter << ','
            << record.chemistryEnergyAbsError << ','
            << record.chemistryEnergyRelError << ','
            << record.velocityMagnitude << ','
            << record.sumYBefore << ','
            << record.sumYAfter << ','
            << record.minYBefore << ','
            << record.minYAfter << ','
            << record.maxClipCorrection << ','
            << record.sumYCorrection << ','
            << record.maxRenormalizationCorrection << ','
            << record.renormalizationFactor << ','
            << record.temperatureClipCorrection << ','
            << record.retrySubsteps << ','
            << record.cvodeInternalSteps << ','
            << record.cvodeReturnFlag << ','
            << record.qdot << ','
            << record.maxAbsWdot << ','
            << record.dominantSpecies[0] << ','
            << record.dominantWdot[0] << ','
            << record.dominantSpecies[1] << ','
            << record.dominantWdot[1] << ','
            << record.dominantSpecies[2] << ','
            << record.dominantWdot[2] << ','
            << '"' << EscapeCSV(SerializeIndexVector(record.clippedSpecies)) << '"'
            << '\n';
      }
    }
  }

  {
    std::ofstream out(samplesPath);
    if (out) {
      out << "rel_error,clip_correction,sumY_correction,renorm_correction\n";
      out << std::scientific << std::setprecision(17);
      for (const auto& sample : ChemEnergyAuditSamples) {
        out << sample.relError << ','
            << sample.clipCorrection << ','
            << sample.sumYCorrection << ','
            << sample.renormalizationCorrection << '\n';
      }
    }
  }
}

unsigned long CDetonationEulerSolver::SetPrimitive_Variables(CConfig* config) {
  unsigned long nonPhysicalPoints = 0;
  for (unsigned long iPoint = 0; iPoint < nPoint; ++iPoint) {
    const bool physical = nodes->SetPrimVar(iPoint, ReactiveModel);
    nodes->SetSecondaryVar(iPoint, ReactiveModel);
    if (!physical) nonPhysicalPoints++;
  }
  return nonPhysicalPoints;
}

void CDetonationEulerSolver::Preprocessing(CGeometry* geometry, CSolver**, CConfig* config,
                                           unsigned short, unsigned short, unsigned short, bool Output) {
  if (Output) return;

  LinSysRes.SetValZero();
  DebugResidualStage("after residual reset", geometry, config);
  ErrorCounter = SetPrimitive_Variables(config);
  DebugPrimitiveStage("preprocessing SetPrimitive", geometry, config);
  if (config->GetComm_Level() == COMM_FULL) {
    unsigned long tmp = ErrorCounter;
    SU2_MPI::Allreduce(&tmp, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
    config->SetNonphysical_Points(ErrorCounter);
  }
}

void CDetonationEulerSolver::SetMax_Eigenvalue(CGeometry* geometry, const CConfig* config) {
  struct SoundSpeed {
    FORCEINLINE su2double operator()(const CDetonationEulerVariable& nodes, unsigned long iPoint,
                                     unsigned long jPoint) const {
      return 0.5 * (nodes.GetSoundSpeed(iPoint) + nodes.GetSoundSpeed(jPoint));
    }
    FORCEINLINE su2double operator()(const CDetonationEulerVariable& nodes, unsigned long iPoint) const {
      return nodes.GetSoundSpeed(iPoint);
    }
  } soundSpeed;

  SetMax_Eigenvalue_impl(soundSpeed, geometry, config);
}

void CDetonationEulerSolver::SetTime_Step(CGeometry* geometry, CSolver** solver_container,
                                          CConfig* config, unsigned short iMesh, unsigned long Iteration) {
  struct SoundSpeed {
    FORCEINLINE su2double operator()(const CDetonationEulerVariable& nodes, unsigned long iPoint,
                                     unsigned long jPoint) const {
      return 0.5 * (nodes.GetSoundSpeed(iPoint) + nodes.GetSoundSpeed(jPoint));
    }
    FORCEINLINE su2double operator()(const CDetonationEulerVariable& nodes, unsigned long iPoint) const {
      return nodes.GetSoundSpeed(iPoint);
    }
  } soundSpeed;

  struct NoViscousLambda {
    FORCEINLINE su2double operator()(const CDetonationEulerVariable&, unsigned long, unsigned long) const { return 0.0; }
    FORCEINLINE su2double operator()(const CDetonationEulerVariable&, unsigned long) const { return 0.0; }
  } lambdaVisc;

  SetTime_Step_impl(soundSpeed, lambdaVisc, geometry, solver_container, config, iMesh, Iteration);
}

void CDetonationEulerSolver::Centered_Residual(CGeometry*, CSolver**, CNumerics**, CConfig*, unsigned short,
                                               unsigned short) {
  SU2_MPI::Error("DETONATION_EULER only supports upwind convective discretization.", CURRENT_FUNCTION);
}

void CDetonationEulerSolver::Upwind_Residual(CGeometry* geometry, CSolver**, CNumerics** numerics,
                                             CConfig* config, unsigned short) {
  CachedNumerics = numerics;
  if (ChemistryOnlyMode(config)) return;
  if (config->GetKind_Upwind_Flow() == UPWIND::DETONATION_HLLC) {
    CUpwHLLC_Detonation::ResetDiagnostics();
  }

  CNumerics* convNumerics = numerics[CONV_TERM];

  for (auto iEdge = 0ul; iEdge < geometry->GetnEdge(); ++iEdge) {
    const auto iPoint = geometry->edges->GetNode(iEdge, 0);
    const auto jPoint = geometry->edges->GetNode(iEdge, 1);

    convNumerics->SetNormal(geometry->edges->GetNormal(iEdge));
    convNumerics->SetPrimitive(nodes->GetPrimitive(iPoint), nodes->GetPrimitive(jPoint));
    convNumerics->SetConservative(nodes->GetSolution(iPoint), nodes->GetSolution(jPoint));

    auto residual = convNumerics->ComputeResidual(config);
    LinSysRes.AddBlock(iPoint, residual);
    LinSysRes.SubtractBlock(jPoint, residual);
  }

  DebugResidualStage("after interior edge residual assembly", geometry, config);
}

void CDetonationEulerSolver::Source_Residual(CGeometry*, CSolver**, CNumerics**, CConfig*, unsigned short) {}

void CDetonationEulerSolver::ApplyExplicitFlowStep(CGeometry* geometry, CSolver** solver_container,
                                                   CConfig* config, su2double dtFactor) {
  if (dtFactor == 1.0) {
    Explicit_Iteration<EULER_EXPLICIT>(geometry, solver_container, config, 0);
    return;
  }

  std::vector<su2double> dtBackup(nPointDomain, 0.0);
  for (unsigned long iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    dtBackup[iPoint] = nodes->GetDelta_Time(iPoint);
    nodes->SetDelta_Time(iPoint, dtFactor * dtBackup[iPoint]);
  }

  Explicit_Iteration<EULER_EXPLICIT>(geometry, solver_container, config, 0);

  for (unsigned long iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    nodes->SetDelta_Time(iPoint, dtBackup[iPoint]);
  }
}

void CDetonationEulerSolver::AssembleFlowResidual(CGeometry* geometry, CSolver** solver_container, CConfig* config) {
  if (CachedNumerics == nullptr) {
    SU2_MPI::Error("DETONATION_EULER STRANG requested flow residual reassembly before numerics were cached.", CURRENT_FUNCTION);
  }

  LinSysRes.SetValZero();
  Upwind_Residual(geometry, solver_container, CachedNumerics, config, MGLevel);

  CNumerics* conv_bound_numerics = CachedNumerics[CONV_BOUND_TERM + omp_get_thread_num() * MAX_TERMS];
  CNumerics* visc_bound_numerics = CachedNumerics[VISC_BOUND_TERM + omp_get_thread_num() * MAX_TERMS];

  for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); ++iMarker) {
    switch (config->GetMarker_All_KindBC(iMarker)) {
      case SUPERSONIC_INLET:
        BC_Supersonic_Inlet(geometry, solver_container, conv_bound_numerics, visc_bound_numerics, config, iMarker);
        break;
      case SUPERSONIC_OUTLET:
        BC_Supersonic_Outlet(geometry, solver_container, conv_bound_numerics, visc_bound_numerics, config, iMarker);
        break;
      case FAR_FIELD:
        BC_Far_Field(geometry, solver_container, conv_bound_numerics, visc_bound_numerics, config, iMarker);
        break;
      default:
        break;
    }
  }

  for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); ++iMarker) {
    switch (config->GetMarker_All_KindBC(iMarker)) {
      case SYMMETRY_PLANE:
        BC_Sym_Plane(geometry, solver_container, conv_bound_numerics, visc_bound_numerics, config, iMarker);
        break;
      case EULER_WALL:
        BC_Euler_Wall(geometry, solver_container, conv_bound_numerics, visc_bound_numerics, config, iMarker);
        break;
      default:
        break;
    }
  }
}

void CDetonationEulerSolver::RecomputeConservativeVector(const su2double* primitive,
                                                         const su2double* massFractions,
                                                         su2double* conservative) const {
  const su2double rho = primitive[nSpecies + nDim + 2];
  std::vector<su2double> y(nSpecies, 0.0);
  for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) y[iSpecies] = massFractions[iSpecies];
  ReactiveModel->NormalizeMassFractions(y, 0.0);
  ReactiveModel->SetMassFractions(y.data(), nSpecies);
  ReactiveModel->SetTDState_rhoT(rho, primitive[nSpecies]);

  for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) conservative[iSpecies] = rho * y[iSpecies];
  su2double velocity2 = 0.0;
  for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
    conservative[nSpecies + iDim] = rho * primitive[nSpecies + 1 + iDim];
    velocity2 += primitive[nSpecies + 1 + iDim] * primitive[nSpecies + 1 + iDim];
  }
  conservative[nVar - 1] = rho * (ReactiveModel->GetStaticEnergy() + 0.5 * velocity2);
}

void CDetonationEulerSolver::BuildBoundaryStateFromPrimitive(const su2double* primitive,
                                                             const su2double* massFractions,
                                                             su2double* conservative) const {
  RecomputeConservativeVector(primitive, massFractions, conservative);
}

void CDetonationEulerSolver::BuildSupersonicInletBoundaryState(CConfig* config, unsigned short,
                                                               std::vector<su2double>& ghostPrimitive,
                                                               std::vector<su2double>& ghostMassFractions,
                                                               std::vector<su2double>& ghostConservative) const {
  for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    ghostMassFractions[iSpecies] = config->GetGas_Composition()[iSpecies];
  }

  const su2double density = config->GetDensity_FreeStreamND();
  for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    ghostPrimitive[iSpecies] = density * ghostMassFractions[iSpecies];
  }
  ghostPrimitive[nSpecies] = config->GetTemperature_FreeStreamND();
  su2double velocity2 = 0.0;
  for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
    ghostPrimitive[nSpecies + 1 + iDim] = config->GetVelocity_FreeStreamND()[iDim];
    velocity2 += ghostPrimitive[nSpecies + 1 + iDim] * ghostPrimitive[nSpecies + 1 + iDim];
  }
  ghostPrimitive[nSpecies + nDim + 1] = config->GetPressure_FreeStreamND();
  ghostPrimitive[nSpecies + nDim + 2] = density;
  ReactiveModel->SetMassFractions(ghostMassFractions.data(), nSpecies);
  ReactiveModel->SetTDState_PT(ghostPrimitive[nSpecies + nDim + 1], ghostPrimitive[nSpecies]);
  ghostPrimitive[nSpecies + nDim + 3] =
      ReactiveModel->GetStaticEnergy() + 0.5 * velocity2 + ghostPrimitive[nSpecies + nDim + 1] / density;
  ghostPrimitive[nSpecies + nDim + 4] = ReactiveModel->GetSoundSpeed();
  BuildBoundaryStateFromPrimitive(ghostPrimitive.data(), ghostMassFractions.data(), ghostConservative.data());
}

void CDetonationEulerSolver::BC_Supersonic_Inlet(CGeometry* geometry, CSolver**, CNumerics* conv_numerics,
                                                 CNumerics*, CConfig* config, unsigned short marker) {
  if (ChemistryOnlyMode(config)) return;

  std::vector<su2double> ghostPrimitive(nPrimVar, 0.0);
  std::vector<su2double> ghostConservative(nVar, 0.0);
  std::vector<su2double> ghostMassFractions(nSpecies, 0.0);
  const bool debugBC = DebugBCResidualEnabled(config);
  BoundaryDebugSummary debugSummary(debugBC, nSpecies, nDim, nVar, nPrimVar);
  BuildSupersonicInletBoundaryState(config, marker, ghostPrimitive, ghostMassFractions, ghostConservative);

  auto applyVertex = [&](unsigned long iVertex) {
    const auto iPoint = geometry->vertex[marker][iVertex]->GetNode();
    if (!geometry->nodes->GetDomain(iPoint)) return;

    su2double normal[MAXNDIM] = {0.0};
    geometry->vertex[marker][iVertex]->GetNormal(normal);
    for (unsigned short iDim = 0; iDim < nDim; ++iDim) normal[iDim] = -normal[iDim];

    conv_numerics->SetNormal(normal);
    conv_numerics->SetPrimitive(nodes->GetPrimitive(iPoint), ghostPrimitive.data());
    conv_numerics->SetConservative(nodes->GetSolution(iPoint), ghostConservative.data());
    auto residual = conv_numerics->ComputeResidual(config);
    LinSysRes.AddBlock(iPoint, residual);
    debugSummary.Accumulate(iVertex, iPoint, geometry->nodes->GetGlobalIndex(iPoint), geometry->nodes->GetCoord(iPoint),
                            normal, nodes->GetPrimitive(iPoint), ghostPrimitive.data(),
                            nodes->GetSolution(iPoint), ghostConservative.data(), residual);
  };

  if (debugBC) {
    for (auto iVertex = 0ul; iVertex < geometry->nVertex[marker]; ++iVertex) applyVertex(iVertex);
  } else {
    SU2_OMP_FOR_DYN(OMP_MIN_SIZE)
    for (auto iVertex = 0ul; iVertex < geometry->nVertex[marker]; ++iVertex) applyVertex(iVertex);
    END_SU2_OMP_FOR
  }

  AppendBCResidualCSV(geometry, config, marker, debugSummary);
  DebugResidualStage("after boundary marker residual assembly", geometry, config,
                     config->GetMarker_All_TagBound(marker));
}

void CDetonationEulerSolver::BC_Supersonic_Outlet(CGeometry* geometry, CSolver**, CNumerics* conv_numerics,
                                                  CNumerics*, CConfig* config, unsigned short marker) {
  if (ChemistryOnlyMode(config)) return;

  std::vector<su2double> ghostConservative(nVar, 0.0);
  const bool debugBC = DebugBCResidualEnabled(config);
  BoundaryDebugSummary debugSummary(debugBC, nSpecies, nDim, nVar, nPrimVar);
  auto applyVertex = [&](unsigned long iVertex) {
    const auto iPoint = geometry->vertex[marker][iVertex]->GetNode();
    if (!geometry->nodes->GetDomain(iPoint)) return;

    su2double normal[MAXNDIM] = {0.0};
    geometry->vertex[marker][iVertex]->GetNormal(normal);
    for (unsigned short iDim = 0; iDim < nDim; ++iDim) normal[iDim] = -normal[iDim];

    for (unsigned short iVar = 0; iVar < nVar; ++iVar) ghostConservative[iVar] = nodes->GetSolution(iPoint, iVar);
    conv_numerics->SetNormal(normal);
    conv_numerics->SetPrimitive(nodes->GetPrimitive(iPoint), nodes->GetPrimitive(iPoint));
    conv_numerics->SetConservative(nodes->GetSolution(iPoint), ghostConservative.data());
    auto residual = conv_numerics->ComputeResidual(config);
    LinSysRes.AddBlock(iPoint, residual);
    debugSummary.Accumulate(iVertex, iPoint, geometry->nodes->GetGlobalIndex(iPoint), geometry->nodes->GetCoord(iPoint),
                            normal, nodes->GetPrimitive(iPoint), nodes->GetPrimitive(iPoint),
                            nodes->GetSolution(iPoint), ghostConservative.data(), residual);
  };

  if (debugBC) {
    for (auto iVertex = 0ul; iVertex < geometry->nVertex[marker]; ++iVertex) applyVertex(iVertex);
  } else {
    SU2_OMP_FOR_DYN(OMP_MIN_SIZE)
    for (auto iVertex = 0ul; iVertex < geometry->nVertex[marker]; ++iVertex) applyVertex(iVertex);
    END_SU2_OMP_FOR
  }

  AppendBCResidualCSV(geometry, config, marker, debugSummary);
  DebugResidualStage("after boundary marker residual assembly", geometry, config,
                     config->GetMarker_All_TagBound(marker));
}

void CDetonationEulerSolver::BC_Far_Field(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                                          CNumerics* visc_numerics, CConfig* config, unsigned short marker) {
  BC_Supersonic_Inlet(geometry, solver_container, conv_numerics, visc_numerics, config, marker);
}

void CDetonationEulerSolver::BC_Sym_Plane(CGeometry* geometry, CSolver**, CNumerics* conv_numerics,
                                          CNumerics*, CConfig* config, unsigned short marker) {
  if (ChemistryOnlyMode(config)) return;

  std::vector<su2double> ghostPrimitive(nPrimVar, 0.0);
  std::vector<su2double> ghostConservative(nVar, 0.0);
  std::vector<su2double> ghostMassFractions(nSpecies, 0.0);
  const bool debugBC = DebugBCResidualEnabled(config);
  BoundaryDebugSummary debugSummary(debugBC, nSpecies, nDim, nVar, nPrimVar);

  auto applyVertex = [&](unsigned long iVertex) {
    const auto iPoint = geometry->vertex[marker][iVertex]->GetNode();
    if (!geometry->nodes->GetDomain(iPoint)) return;

    for (unsigned short iVar = 0; iVar < nPrimVar; ++iVar) ghostPrimitive[iVar] = nodes->GetPrimitive(iPoint, iVar);
    for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) ghostMassFractions[iSpecies] = nodes->GetMassFraction(iPoint, iSpecies);

    su2double normal[MAXNDIM] = {0.0};
    su2double unitNormal[MAXNDIM] = {0.0};
    geometry->vertex[marker][iVertex]->GetNormal(normal);
    const su2double area = GeometryToolbox::Norm(nDim, normal);
    for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
      normal[iDim] = -normal[iDim];
      unitNormal[iDim] = normal[iDim] / area;
    }

    su2double vn = 0.0;
    for (unsigned short iDim = 0; iDim < nDim; ++iDim) vn += ghostPrimitive[nSpecies + 1 + iDim] * unitNormal[iDim];
    for (unsigned short iDim = 0; iDim < nDim; ++iDim) ghostPrimitive[nSpecies + 1 + iDim] -= 2.0 * vn * unitNormal[iDim];
    BuildBoundaryStateFromPrimitive(ghostPrimitive.data(), ghostMassFractions.data(), ghostConservative.data());

    conv_numerics->SetNormal(normal);
    conv_numerics->SetPrimitive(nodes->GetPrimitive(iPoint), ghostPrimitive.data());
    conv_numerics->SetConservative(nodes->GetSolution(iPoint), ghostConservative.data());
    auto residual = conv_numerics->ComputeResidual(config);
    LinSysRes.AddBlock(iPoint, residual);
    debugSummary.Accumulate(iVertex, iPoint, geometry->nodes->GetGlobalIndex(iPoint), geometry->nodes->GetCoord(iPoint),
                            normal, nodes->GetPrimitive(iPoint), ghostPrimitive.data(),
                            nodes->GetSolution(iPoint), ghostConservative.data(), residual);
  };

  if (debugBC) {
    for (auto iVertex = 0ul; iVertex < geometry->nVertex[marker]; ++iVertex) applyVertex(iVertex);
  } else {
    SU2_OMP_FOR_DYN(OMP_MIN_SIZE)
    for (auto iVertex = 0ul; iVertex < geometry->nVertex[marker]; ++iVertex) applyVertex(iVertex);
    END_SU2_OMP_FOR
  }

  AppendBCResidualCSV(geometry, config, marker, debugSummary);
  DebugResidualStage("after boundary marker residual assembly", geometry, config,
                     config->GetMarker_All_TagBound(marker));
}

void CDetonationEulerSolver::AdvanceChemistry(CGeometry* geometry, CConfig* config) {
  if (config->GetFrozen()) {
    ChemFailCells = 0;
    ChemSubstepsMax = 0;
    ChemEnergyAbsMax = 0.0;
    ChemEnergyRelMax = 0.0;
    std::fill(ChemQdot.begin(), ChemQdot.end(), 0.0);
    std::fill(ChemFailFlag.begin(), ChemFailFlag.end(), 0u);
    std::fill(ChemSubstepsLast.begin(), ChemSubstepsLast.end(), 0ul);
    return;
  }

  ChemFailCells = 0;
  ChemSubstepsMax = 0;
  ChemEnergyAbsMax = 0.0;
  ChemEnergyRelMax = 0.0;
  std::fill(ChemQdot.begin(), ChemQdot.end(), 0.0);
  std::fill(ChemFailFlag.begin(), ChemFailFlag.end(), 0u);
  std::fill(ChemSubstepsLast.begin(), ChemSubstepsLast.end(), 0ul);

  const bool auditEnabled = ChemistryEnergyAuditEnabled(config);
  bool haveLocalWorst = false;
  ChemEnergyWorstRecord localWorstRecord;
  localWorstRecord.chemistryEnergyRelError = -1.0;

  for (unsigned long iPoint = 0; iPoint < nPointDomain; ++iPoint) {
    const su2double rho = nodes->GetDensity(iPoint);
    std::vector<su2double> y(nSpecies, 0.0);
    for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) y[iSpecies] = nodes->GetMassFraction(iPoint, iSpecies);

    const auto result = ChemistryIntegrator->IntegrateCell(nodes->GetDelta_Time(iPoint), rho, nodes->GetTemperature(iPoint), y);
    ChemSubstepsLast[iPoint] = result.substeps;
    if (!result.success) {
      ChemFailCells++;
      ChemFailFlag[iPoint] = 1u;
      ChemFailCount[iPoint] += 1ul;
      continue;
    }

    ChemSubstepsMax = max(ChemSubstepsMax, result.substeps);
    ChemEnergyAbsMax = max(ChemEnergyAbsMax, fabs(result.chemistryEnergyAbsError));
    ChemEnergyRelMax = max(ChemEnergyRelMax, result.chemistryEnergyRelError);
    ChemQdot[iPoint] = result.qdot;

    if (auditEnabled) {
      ChemEnergyAuditSamples.push_back(
          {result.chemistryEnergyRelError, result.maxClipCorrection, result.sumYCorrection,
           result.maxRenormalizationCorrection});

      if (!haveLocalWorst || result.chemistryEnergyRelError > localWorstRecord.chemistryEnergyRelError) {
        haveLocalWorst = true;
        localWorstRecord.timeIter = config->GetTimeIter();
        localWorstRecord.physicalTime = config->GetPhysicalTime();
        localWorstRecord.globalNode = geometry->nodes->GetGlobalIndex(iPoint);
        const auto coord = geometry->nodes->GetCoord(iPoint);
        for (unsigned short iDim = 0; iDim < nDim; ++iDim) localWorstRecord.coord[iDim] = coord[iDim];
        localWorstRecord.density = rho;
        localWorstRecord.temperatureBefore = result.temperatureBefore;
        localWorstRecord.temperatureAfter = result.temperature;
        localWorstRecord.pressureBefore = result.pressureBefore;
        localWorstRecord.pressureAfter = result.pressureAfter;
        localWorstRecord.internalEnergyBefore = result.internalEnergyBefore;
        localWorstRecord.internalEnergyAfter = result.internalEnergyAfter;
        localWorstRecord.chemistryEnergyAbsError = result.chemistryEnergyAbsError;
        localWorstRecord.chemistryEnergyRelError = result.chemistryEnergyRelError;
        su2double velocity2 = 0.0;
        for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
          const su2double vel = nodes->GetVelocity(iPoint, iDim);
          velocity2 += vel * vel;
        }
        localWorstRecord.velocityMagnitude = sqrt(velocity2);
        localWorstRecord.sumYBefore = result.sumYBefore;
        localWorstRecord.sumYAfter = result.sumYAfter;
        localWorstRecord.minYBefore = result.minYBefore;
        localWorstRecord.minYAfter = result.minYAfter;
        localWorstRecord.maxClipCorrection = result.maxClipCorrection;
        localWorstRecord.sumYCorrection = result.sumYCorrection;
        localWorstRecord.maxRenormalizationCorrection = result.maxRenormalizationCorrection;
        localWorstRecord.renormalizationFactor = result.renormalizationFactor;
        localWorstRecord.temperatureClipCorrection = result.temperatureClipCorrection;
        localWorstRecord.retrySubsteps = result.substeps;
        localWorstRecord.cvodeInternalSteps = result.internalSteps;
        localWorstRecord.cvodeReturnFlag = result.cvodeReturnFlag;
        localWorstRecord.qdot = result.qdot;
        localWorstRecord.maxAbsWdot = result.maxAbsWdot;
        localWorstRecord.dominantSpecies = result.dominantSpecies;
        localWorstRecord.dominantWdot = result.dominantWdot;
        localWorstRecord.clippedSpecies = result.clippedSpecies;
      }
    }

    std::vector<su2double> primitive(nPrimVar, 0.0);
    for (unsigned short iVar = 0; iVar < nPrimVar; ++iVar) primitive[iVar] = nodes->GetPrimitive(iPoint, iVar);
    primitive[nSpecies] = result.temperature;
    std::vector<su2double> conservative(nVar, 0.0);
    RecomputeConservativeVector(primitive.data(), result.massFractions.data(), conservative.data());
    for (unsigned short iVar = 0; iVar < nVar; ++iVar) nodes->SetSolution(iPoint, iVar, conservative[iVar]);
  }

  if (config->GetComm_Level() == COMM_FULL) {
    unsigned long failCellsGlobal = ChemFailCells;
    unsigned long substepsMaxGlobal = ChemSubstepsMax;
    su2double energyAbsMaxGlobal = ChemEnergyAbsMax;
    su2double energyRelMaxGlobal = ChemEnergyRelMax;

    SU2_MPI::Allreduce(&failCellsGlobal, &ChemFailCells, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
    SU2_MPI::Allreduce(&substepsMaxGlobal, &ChemSubstepsMax, 1, MPI_UNSIGNED_LONG, MPI_MAX, SU2_MPI::GetComm());
    SU2_MPI::Allreduce(&energyAbsMaxGlobal, &ChemEnergyAbsMax, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());
    SU2_MPI::Allreduce(&energyRelMaxGlobal, &ChemEnergyRelMax, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());
  }

  if (auditEnabled && haveLocalWorst) {
    ChemEnergyWorstRecords.push_back(localWorstRecord);
  }

  if (SU2_MPI::GetRank() == MASTER_NODE && ChemEnergyRelMax > 1.0e-8) {
    cout << "WARNING: DETONATION_EULER chemistry energy closure drift exceeded 1e-8 in this step. "
         << "CHEM_E_REL_MAX = " << ChemEnergyRelMax
         << " (treat verification as failed if this exceeds 1e-6)." << endl;
  }

  InitiateComms(geometry, config, MPI_QUANTITIES::SOLUTION);
  CompleteComms(geometry, config, MPI_QUANTITIES::SOLUTION);
}

void CDetonationEulerSolver::ExplicitEuler_Iteration(CGeometry* geometry, CSolver** solver_container, CConfig* config) {
  SaveDebugSolutionSnapshot();
  DebugResidualStage("after MPI residual communication", geometry, config);
  DebugResidualStage("before explicit update", geometry, config);

  const bool useStrang = UseStrangSplitting(config) && !ChemistryOnlyMode(config) && !config->GetFrozen();

  if (useStrang) {
    nodes->Set_OldSolution();
    ApplyExplicitFlowStep(geometry, solver_container, config, 0.5);
  } else if (!ChemistryOnlyMode(config)) {
    ApplyExplicitFlowStep(geometry, solver_container, config, 1.0);
  } else if (DebugResidualBreakdownEnabled(config) && SU2_MPI::GetRank() == MASTER_NODE) {
    cout << "WARNING: DETONATION_CHEMISTRY_ONLY is for 0D verification only." << endl;
  }

  DebugResidualStage("after explicit update", geometry, config);
  SetPrimitive_Variables(config);
  DebugResidualStage("after SetPrimitive before chemistry", geometry, config);
  DebugPrimitiveStage("after SetPrimitive before chemistry", geometry, config);
  AdvanceChemistry(geometry, config);
  DebugResidualStage("after chemistry", geometry, config);
  SetPrimitive_Variables(config);
  DebugResidualStage("after SetPrimitive after chemistry", geometry, config);
  DebugPrimitiveStage("after SetPrimitive after chemistry", geometry, config);

  if (useStrang) {
    nodes->Set_OldSolution();
    SaveDebugSolutionSnapshot();
    AssembleFlowResidual(geometry, solver_container, config);
    DebugResidualStage("after Strang second flow residual assembly", geometry, config);
    ApplyExplicitFlowStep(geometry, solver_container, config, 0.5);
    DebugResidualStage("after Strang second flow half-step", geometry, config);
    SetPrimitive_Variables(config);
    DebugResidualStage("after SetPrimitive after Strang", geometry, config);
    DebugPrimitiveStage("after SetPrimitive after Strang", geometry, config);
  }
}
