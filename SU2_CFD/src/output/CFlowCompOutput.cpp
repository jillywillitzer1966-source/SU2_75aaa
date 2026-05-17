/*!
 * \file CFlowCompOutput.cpp
 * \brief Main subroutines for compressible flow output
 * \author R. Sanchez
 * \version 8.1.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2024, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../../include/output/CFlowCompOutput.hpp"

#include <limits>
#include <vector>

#include "../../../Common/include/geometry/CGeometry.hpp"
#include "../../include/solvers/CSolver.hpp"
#include "../../include/solvers/CDetonationEulerSolver.hpp"
#include "../../include/variables/CDetonationEulerVariable.hpp"

CFlowCompOutput::CFlowCompOutput(const CConfig *config, unsigned short nDim) : CFlowOutput(config, nDim, false) {
  const bool detonation = config->GetKind_Solver() == MAIN_SOLVER::DETONATION_EULER;

  turb_model = config->GetKind_Turb_Model();

  /*--- Set the default history fields if nothing is set in the config file ---*/

  if (nRequestedHistoryFields == 0){
    RequestCommonHistory(config->GetTime_Domain());
    if (detonation) {
      requestedHistoryFields.emplace_back("CHEMISTRY");
      requestedHistoryFields.emplace_back("DETONATION_HLLC");
    }
    nRequestedHistoryFields = requestedHistoryFields.size();
  }
  if (nRequestedScreenFields == 0){
    if (config->GetTime_Domain()) requestedScreenFields.emplace_back("TIME_ITER");
    if (multiZone) requestedScreenFields.emplace_back("OUTER_ITER");
    requestedScreenFields.emplace_back("INNER_ITER");
    requestedScreenFields.emplace_back("RMS_DENSITY");
    requestedScreenFields.emplace_back("RMS_MOMENTUM-X");
    requestedScreenFields.emplace_back("RMS_MOMENTUM-Y");
    requestedScreenFields.emplace_back("RMS_ENERGY");
    if (detonation) {
      requestedScreenFields.emplace_back("CHEM_FAIL_CELLS");
      requestedScreenFields.emplace_back("CHEM_SUBSTEPS_MAX");
    }
    nRequestedScreenFields = requestedScreenFields.size();
  }
  if (nRequestedVolumeFields == 0){
    requestedVolumeFields.emplace_back("COORDINATES");
    requestedVolumeFields.emplace_back("SOLUTION");
    requestedVolumeFields.emplace_back("PRIMITIVE");
    if (detonation) requestedVolumeFields.emplace_back("AUXILIARY");
    nRequestedVolumeFields = requestedVolumeFields.size();
  }

  if (gridMovement) {
    auto notFound = requestedVolumeFields.end();
    if (find(requestedVolumeFields.begin(), notFound, string("GRID_VELOCITY")) == notFound) {
      requestedVolumeFields.emplace_back("GRID_VELOCITY");
      nRequestedVolumeFields++;
    }
  }

  stringstream ss;
  ss << "Zone " << config->GetiZone() << " (Comp. Fluid)";
  multiZoneHeaderString = ss.str();

  /*--- Set the volume filename --- */

  volumeFilename = config->GetVolume_FileName();

  /*--- Set the surface filename --- */

  surfaceFilename = config->GetSurfCoeff_FileName();

  /*--- Set the restart filename --- */

  restartFilename = config->GetRestart_FileName();

  /*--- Set the default convergence field --- */

  if (convFields.empty() ) convFields.emplace_back("RMS_DENSITY");

  if (config->GetFixed_CL_Mode()) {
    if (std::find(convFields.begin(), convFields.end(), "LIFT") != convFields.end()) {
      if (rank == MASTER_NODE)
        cout<<"  Fixed CL: Adding LIFT as Convergence Field to ensure convergence to target CL"<<endl;
      convFields.emplace_back("LIFT");
      newFunc.resize(convFields.size());
      oldFunc.resize(convFields.size());
      cauchySerie.resize(convFields.size(), vector<su2double>(nCauchy_Elems, 0.0));
    }
  }
}

void CFlowCompOutput::SetHistoryOutputFields(CConfig *config){

  /// BEGIN_GROUP: RMS_RES, DESCRIPTION: The root-mean-square residuals of the SOLUTION variables.
  /// DESCRIPTION: Root-mean square residual of the density.
  AddHistoryOutput("RMS_DENSITY",    "rms[Rho]",  ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the density.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the momentum x-component.
  AddHistoryOutput("RMS_MOMENTUM-X", "rms[RhoU]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the momentum x-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the momentum y-component.
  AddHistoryOutput("RMS_MOMENTUM-Y", "rms[RhoV]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the momentum y-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the momentum z-component.
  if (nDim == 3) AddHistoryOutput("RMS_MOMENTUM-Z", "rms[RhoW]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the momentum z-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Root-mean square residual of the energy.
  AddHistoryOutput("RMS_ENERGY",     "rms[RhoE]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the energy.", HistoryFieldType::RESIDUAL);

  AddHistoryOutputFields_ScalarRMS_RES(config);
  /// END_GROUP

  /// BEGIN_GROUP: MAX_RES, DESCRIPTION: The maximum residuals of the SOLUTION variables.
  /// DESCRIPTION: Maximum residual of the density.
  AddHistoryOutput("MAX_DENSITY",    "max[Rho]",  ScreenOutputFormat::FIXED,   "MAX_RES", "Maximum square residual of the density.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the momentum x-component.
  AddHistoryOutput("MAX_MOMENTUM-X", "max[RhoU]", ScreenOutputFormat::FIXED,   "MAX_RES", "Maximum square residual of the momentum x-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the momentum y-component.
  AddHistoryOutput("MAX_MOMENTUM-Y", "max[RhoV]", ScreenOutputFormat::FIXED,   "MAX_RES", "Maximum square residual of the momentum y-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the momentum z-component.
  if (nDim == 3) AddHistoryOutput("MAX_MOMENTUM-Z", "max[RhoW]", ScreenOutputFormat::FIXED,"MAX_RES", "Maximum residual of the z-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the energy.
  AddHistoryOutput("MAX_ENERGY",     "max[RhoE]", ScreenOutputFormat::FIXED,   "MAX_RES", "Maximum residual of the energy.", HistoryFieldType::RESIDUAL);

  AddHistoryOutputFields_ScalarMAX_RES(config);
  /// END_GROUP

  /// BEGIN_GROUP: BGS_RES, DESCRIPTION: The block Gauss Seidel residuals of the SOLUTION variables.
  /// DESCRIPTION: Maximum residual of the density.
  AddHistoryOutput("BGS_DENSITY",    "bgs[Rho]",  ScreenOutputFormat::FIXED,   "BGS_RES", "BGS residual of the density.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the momentum x-component.
  AddHistoryOutput("BGS_MOMENTUM-X", "bgs[RhoU]", ScreenOutputFormat::FIXED,   "BGS_RES", "BGS residual of the momentum x-component.", HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the momentum y-component.
  AddHistoryOutput("BGS_MOMENTUM-Y", "bgs[RhoV]", ScreenOutputFormat::FIXED,   "BGS_RES", "BGS residual of the momentum y-component.",  HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the momentum z-component.
  if (nDim == 3) AddHistoryOutput("BGS_MOMENTUM-Z", "bgs[RhoW]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the z-component.",  HistoryFieldType::RESIDUAL);
  /// DESCRIPTION: Maximum residual of the energy.
  AddHistoryOutput("BGS_ENERGY",     "bgs[RhoE]", ScreenOutputFormat::FIXED,   "BGS_RES", "BGS residual of the energy.",  HistoryFieldType::RESIDUAL);

  AddHistoryOutputFields_ScalarBGS_RES(config);
  /// END_GROUP

  vector<string> Marker_Monitoring;
  for (unsigned short iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++){
    Marker_Monitoring.push_back(config->GetMarker_Monitoring_TagBound(iMarker_Monitoring));
  }

  if (config->GetAeroelastic_Simulation()) {
    /// BEGIN_GROUP: AEROELASTIC, DESCRIPTION: Aeroelastic plunge, pitch
    /// DESCRIPTION: Aeroelastic plunge
    AddHistoryOutputPerSurface("PLUNGE", "plunge", ScreenOutputFormat::FIXED, "AEROELASTIC", Marker_Monitoring, HistoryFieldType::COEFFICIENT);
    /// DESCRIPTION: Aeroelastic pitch
    AddHistoryOutputPerSurface("PITCH", "pitch", ScreenOutputFormat::FIXED, "AEROELASTIC", Marker_Monitoring, HistoryFieldType::COEFFICIENT);
    /// END_GROUP
  }

  /// DESCRIPTION: Linear solver iterations
  AddHistoryOutput("LINSOL_ITER", "Linear_Solver_Iterations", ScreenOutputFormat::INTEGER, "LINSOL", "Number of iterations of the linear solver.");
  AddHistoryOutput("LINSOL_RESIDUAL", "LinSolRes", ScreenOutputFormat::FIXED, "LINSOL", "Residual of the linear solver.");
  AddHistoryOutputFieldsScalarLinsol(config);

  AddHistoryOutput("MIN_DELTA_TIME", "Min DT", ScreenOutputFormat::SCIENTIFIC, "CFL_NUMBER", "Current minimum local time step");
  AddHistoryOutput("MAX_DELTA_TIME", "Max DT", ScreenOutputFormat::SCIENTIFIC, "CFL_NUMBER", "Current maximum local time step");

  AddHistoryOutput("MIN_CFL", "Min CFL", ScreenOutputFormat::SCIENTIFIC, "CFL_NUMBER", "Current minimum of the local CFL numbers");
  AddHistoryOutput("MAX_CFL", "Max CFL", ScreenOutputFormat::SCIENTIFIC, "CFL_NUMBER", "Current maximum of the local CFL numbers");
  AddHistoryOutput("AVG_CFL", "Avg CFL", ScreenOutputFormat::SCIENTIFIC, "CFL_NUMBER", "Current average of the local CFL numbers");

  if (config->GetKind_Solver() == MAIN_SOLVER::DETONATION_EULER) {
    AddHistoryOutput("CHEM_FAIL_CELLS", "ChemFailCells", ScreenOutputFormat::INTEGER, "CHEMISTRY",
                     "Number of cells whose chemistry substep failed in the current outer iteration.");
    AddHistoryOutput("CHEM_SUBSTEPS_MAX", "ChemSubstepsMax", ScreenOutputFormat::INTEGER, "CHEMISTRY",
                     "Maximum number of chemistry retry substeps used by any cell in the current outer iteration.");
    AddHistoryOutput("CHEM_E_ABS_MAX", "ChemEAbsMax", ScreenOutputFormat::SCIENTIFIC, "CHEMISTRY",
                     "Maximum absolute internal-energy drift across constant-volume chemistry substeps.");
    AddHistoryOutput("CHEM_E_REL_MAX", "ChemERelMax", ScreenOutputFormat::SCIENTIFIC, "CHEMISTRY",
                     "Maximum relative internal-energy drift across constant-volume chemistry substeps.");
    AddHistoryOutput("HLLC_FALLBACK_COUNT", "HLLCFallback", ScreenOutputFormat::INTEGER, "DETONATION_HLLC",
                     "Number of DETONATION_HLLC interfaces that fell back to DETONATION_LLF in the current residual assembly.");
    AddHistoryOutput("HLLC_NEGATIVE_STATE_COUNT", "HLLCNegativeState", ScreenOutputFormat::INTEGER, "DETONATION_HLLC",
                     "Number of DETONATION_HLLC interfaces with rejected non-positive candidate states.");
    AddHistoryOutput("HLLC_MIN_PRESSURE", "HLLCMinPressure", ScreenOutputFormat::SCIENTIFIC, "DETONATION_HLLC",
                     "Minimum pressure sampled by DETONATION_HLLC diagnostics in the current residual assembly.");
    AddHistoryOutput("HLLC_MIN_DENSITY", "HLLCMinDensity", ScreenOutputFormat::SCIENTIFIC, "DETONATION_HLLC",
                     "Minimum density sampled by DETONATION_HLLC diagnostics in the current residual assembly.");
    AddHistoryOutput("DET_AVG_TEMPERATURE", "DetAvgTemperature", ScreenOutputFormat::SCIENTIFIC, "DETONATION_STATE",
                     "Domain-averaged temperature for DETONATION_EULER.");
    AddHistoryOutput("DET_AVG_PRESSURE", "DetAvgPressure", ScreenOutputFormat::SCIENTIFIC, "DETONATION_STATE",
                     "Domain-averaged pressure for DETONATION_EULER.");
    AddHistoryOutput("DET_AVG_DENSITY", "DetAvgDensity", ScreenOutputFormat::SCIENTIFIC, "DETONATION_STATE",
                     "Domain-averaged density for DETONATION_EULER.");
    AddHistoryOutput("DET_AVG_INTERNAL_ENERGY", "DetAvgInternalEnergy", ScreenOutputFormat::SCIENTIFIC, "DETONATION_STATE",
                     "Domain-averaged specific internal energy for DETONATION_EULER.");
    AddHistoryOutput("DET_AVG_QDOT", "DetAvgQdot", ScreenOutputFormat::SCIENTIFIC, "DETONATION_STATE",
                     "Domain-averaged diagnostic heat-release proxy for DETONATION_EULER.");
    AddHistoryOutput("DET_SPAN_TEMPERATURE", "DetSpanTemperature", ScreenOutputFormat::SCIENTIFIC, "HOMOGENEITY",
                     "Spatial max-min temperature spread for DETONATION_EULER.");
    AddHistoryOutput("DET_SPAN_PRESSURE", "DetSpanPressure", ScreenOutputFormat::SCIENTIFIC, "HOMOGENEITY",
                     "Spatial max-min pressure spread for DETONATION_EULER.");
    AddHistoryOutput("DET_SPAN_DENSITY", "DetSpanDensity", ScreenOutputFormat::SCIENTIFIC, "HOMOGENEITY",
                     "Spatial max-min density spread for DETONATION_EULER.");
    AddHistoryOutput("DET_SPAN_MOMENTUM-X", "DetSpanMomentumX", ScreenOutputFormat::SCIENTIFIC, "HOMOGENEITY",
                     "Spatial max-min x-momentum spread for DETONATION_EULER.");
    AddHistoryOutput("DET_SPAN_MOMENTUM-Y", "DetSpanMomentumY", ScreenOutputFormat::SCIENTIFIC, "HOMOGENEITY",
                     "Spatial max-min y-momentum spread for DETONATION_EULER.");
    if (nDim == 3) {
      AddHistoryOutput("DET_SPAN_MOMENTUM-Z", "DetSpanMomentumZ", ScreenOutputFormat::SCIENTIFIC, "HOMOGENEITY",
                       "Spatial max-min z-momentum spread for DETONATION_EULER.");
    }
    AddHistoryOutput("DET_SPAN_ENERGY", "DetSpanEnergy", ScreenOutputFormat::SCIENTIFIC, "HOMOGENEITY",
                     "Spatial max-min total-energy spread for DETONATION_EULER.");
    AddHistoryOutput("DET_MIN_MASS_FRACTION", "DetMinMassFraction", ScreenOutputFormat::SCIENTIFIC, "HOMOGENEITY",
                     "Minimum species mass fraction across all cells and species.");
    AddHistoryOutput("DET_MAX_SUM_Y_ERROR", "DetMaxSumYError", ScreenOutputFormat::SCIENTIFIC, "HOMOGENEITY",
                     "Maximum absolute deviation of sum(Y) from unity across all cells.");
    for (unsigned short iSpecies = 0; iSpecies < config->GetnSpecies(); ++iSpecies) {
      const auto suffix = std::to_string(iSpecies);
      AddHistoryOutput("DET_AVG_Y_" + suffix, "DetAvgY_" + suffix, ScreenOutputFormat::SCIENTIFIC, "DETONATION_STATE",
                       "Domain-averaged species mass fraction " + suffix + " for DETONATION_EULER.");
      AddHistoryOutput("DET_SPAN_Y_" + suffix, "DetSpanY_" + suffix, ScreenOutputFormat::SCIENTIFIC, "HOMOGENEITY",
                       "Spatial max-min species mass fraction spread for species " + suffix + '.');
      AddHistoryOutput("DET_SPAN_RHOS_" + suffix, "DetSpanRhoS_" + suffix, ScreenOutputFormat::SCIENTIFIC, "HOMOGENEITY",
                       "Spatial max-min species partial-density spread for species " + suffix + '.');
    }
  }

  /// BEGIN_GROUP: FIXED_CL, DESCRIPTION: Relevant outputs for the Fixed CL mode
  if (config->GetFixed_CL_Mode()){
    /// DESCRIPTION: Difference between current and target CL
    AddHistoryOutput("DELTA_CL", "Delta_CL", ScreenOutputFormat::SCIENTIFIC, "FIXED_CL", "Difference between Target CL and current CL", HistoryFieldType::COEFFICIENT);
    /// DESCRIPTION: Angle of attack before the most recent update
    AddHistoryOutput("PREV_AOA", "Previous_AOA", ScreenOutputFormat::FIXED, "FIXED_CL", "Angle of Attack at the previous iteration of the Fixed CL driver");
    /// DESCRIPTION: Last change in angle of attack by the Fixed CL driver
    AddHistoryOutput("CHANGE_IN_AOA", "Change_in_AOA", ScreenOutputFormat::SCIENTIFIC, "FIXED_CL", "Last change in Angle of Attack by Fixed CL Driver", HistoryFieldType::RESIDUAL);
    /// DESCRIPTION: AOA control command by the CL Driver
    AddHistoryOutput("CL_DRIVER_COMMAND", "CL_Driver_Command", ScreenOutputFormat::SCIENTIFIC, "FIXED_CL", "CL Driver's control command", HistoryFieldType::RESIDUAL);
  }
  /// END_GROUP

  if (config->GetDeform_Mesh()){
    AddHistoryOutput("DEFORM_MIN_VOLUME", "MinVolume", ScreenOutputFormat::SCIENTIFIC, "DEFORM", "Minimum volume in the mesh");
    AddHistoryOutput("DEFORM_MAX_VOLUME", "MaxVolume", ScreenOutputFormat::SCIENTIFIC, "DEFORM", "Maximum volume in the mesh");
    AddHistoryOutput("DEFORM_ITER", "DeformIter", ScreenOutputFormat::INTEGER, "DEFORM", "Linear solver iterations for the mesh deformation");
    AddHistoryOutput("DEFORM_RESIDUAL", "DeformRes", ScreenOutputFormat::FIXED, "DEFORM", "Residual of the linear solver for the mesh deformation");
  }

  AddAnalyzeSurfaceOutput(config);

  AddAerodynamicCoefficients(config);

  if (config->GetViscous()) {
    AddHistoryOutput("BUFFET", "Buffet", ScreenOutputFormat::SCIENTIFIC, "AERO_COEFF", "Buffet sensor", HistoryFieldType::COEFFICIENT);
  }

  AddHeatCoefficients(config);

  AddRotatingFrameCoefficients();

  AddCpInverseDesignOutput();

  AddNearfieldInverseDesignOutput();

  if (config->GetBoolTurbomachinery()) AddTurboOutput(config->GetnZone());

}

void CFlowCompOutput::SetVolumeOutputFields(CConfig *config){
  const bool detonation = config->GetKind_Solver() == MAIN_SOLVER::DETONATION_EULER;
  const auto nSpecies = config->GetnSpecies();

  // Grid coordinates
  AddCoordinates();

  // Solution variables
  AddVolumeOutput("DENSITY",    "Density",    "SOLUTION", "Density");
  AddVolumeOutput("MOMENTUM-X", "Momentum_x", "SOLUTION", "x-component of the momentum vector");
  AddVolumeOutput("MOMENTUM-Y", "Momentum_y", "SOLUTION", "y-component of the momentum vector");

  if (nDim == 3)
    AddVolumeOutput("MOMENTUM-Z", "Momentum_z", "SOLUTION", "z-component of the momentum vector");
  AddVolumeOutput("ENERGY",     "Energy",     "SOLUTION", "Energy");

  SetVolumeOutputFieldsScalarSolution(config);

  if (detonation) {
    for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      AddVolumeOutput("RHOS_" + std::to_string(iSpecies), "RhoS_" + std::to_string(iSpecies), "AUXILIARY",
                      "Species partial density " + std::to_string(iSpecies));
      AddVolumeOutput("MASSFRAC_" + std::to_string(iSpecies), "MassFrac_" + std::to_string(iSpecies), "AUXILIARY",
                      "Species mass fraction " + std::to_string(iSpecies));
    }
    AddVolumeOutput("QDOT", "Qdot", "AUXILIARY", "Volumetric chemical heat release source proxy.");
    AddVolumeOutput("CHEM_FAIL_FLAG", "Chem_Fail_Flag", "AUXILIARY",
                    "Per-step chemistry failure flag for the current cell.");
    AddVolumeOutput("CHEM_FAIL_COUNT", "Chem_Fail_Count", "AUXILIARY",
                    "Cumulative number of chemistry failures for the current cell.");
    AddVolumeOutput("CHEM_SUBSTEPS", "Chem_Substeps", "AUXILIARY",
                    "Number of chemistry substeps used in the most recent chemistry advance for the current cell.");
  }

  // Grid velocity
  if (gridMovement){
    AddVolumeOutput("GRID_VELOCITY-X", "Grid_Velocity_x", "GRID_VELOCITY", "x-component of the grid velocity vector");
    AddVolumeOutput("GRID_VELOCITY-Y", "Grid_Velocity_y", "GRID_VELOCITY", "y-component of the grid velocity vector");
    if (nDim == 3 )
      AddVolumeOutput("GRID_VELOCITY-Z", "Grid_Velocity_z", "GRID_VELOCITY", "z-component of the grid velocity vector");
  }

  // Primitive variables
  AddVolumeOutput("PRESSURE",    "Pressure",                "PRIMITIVE", "Pressure");
  AddVolumeOutput("TEMPERATURE", "Temperature",             "PRIMITIVE", "Temperature");
  AddVolumeOutput("MACH",        "Mach",                    "PRIMITIVE", "Mach number");
  AddVolumeOutput("SOUND_SPEED", "Sound_Speed",             "PRIMITIVE", "Frozen-state sound speed");
  AddVolumeOutput("PRESSURE_COEFF", "Pressure_Coefficient", "PRIMITIVE", "Pressure coefficient");
  AddVolumeOutput("VELOCITY-X", "Velocity_x", "PRIMITIVE", "x-component of the velocity vector");
  AddVolumeOutput("VELOCITY-Y", "Velocity_y", "PRIMITIVE", "y-component of the velocity vector");

  if (nDim == 3)
    AddVolumeOutput("VELOCITY-Z", "Velocity_z", "PRIMITIVE", "z-component of the velocity vector");

  // Datadriven fluid model
  if(config->GetKind_FluidModel() == DATADRIVEN_FLUID){
    AddVolumeOutput("EXTRAPOLATION", "Extrapolation", "PRIMITIVE", "Density, energy outside data range");
    AddVolumeOutput("FLUIDMODEL_NEWTONITER", "nIter_Newton", "PRIMITIVE", "Number of iterations evaluated by the Newton solver");
    AddVolumeOutput("ENTROPY", "Entropy", "PRIMITIVE", "Fluid entropy value");
  }

  if (config->GetViscous()) {
    AddVolumeOutput("LAMINAR_VISCOSITY", "Laminar_Viscosity", "PRIMITIVE", "Laminar viscosity");

    AddVolumeOutput("SKIN_FRICTION-X", "Skin_Friction_Coefficient_x", "PRIMITIVE", "x-component of the skin friction vector");
    AddVolumeOutput("SKIN_FRICTION-Y", "Skin_Friction_Coefficient_y", "PRIMITIVE", "y-component of the skin friction vector");
    if (nDim == 3)
      AddVolumeOutput("SKIN_FRICTION-Z", "Skin_Friction_Coefficient_z", "PRIMITIVE", "z-component of the skin friction vector");

    AddVolumeOutput("HEAT_FLUX", "Heat_Flux", "PRIMITIVE", "Heat-flux");
    AddVolumeOutput("Y_PLUS", "Y_Plus", "PRIMITIVE", "Non-dim. wall distance (Y-Plus)");
  }

  SetVolumeOutputFieldsScalarPrimitive(config);

  //Residuals
  AddVolumeOutput("RES_DENSITY", "Residual_Density", "RESIDUAL", "Residual of the density");
  AddVolumeOutput("RES_MOMENTUM-X", "Residual_Momentum_x", "RESIDUAL", "Residual of the x-momentum component");
  AddVolumeOutput("RES_MOMENTUM-Y", "Residual_Momentum_y", "RESIDUAL", "Residual of the y-momentum component");
  if (nDim == 3)
    AddVolumeOutput("RES_MOMENTUM-Z", "Residual_Momentum_z", "RESIDUAL", "Residual of the z-momentum component");
  AddVolumeOutput("RES_ENERGY", "Residual_Energy", "RESIDUAL", "Residual of the energy");

  SetVolumeOutputFieldsScalarResidual(config);

  if (config->GetKind_SlopeLimit_Flow() != LIMITER::NONE && config->GetKind_SlopeLimit_Flow() != LIMITER::VAN_ALBADA_EDGE) {
    AddVolumeOutput("LIMITER_VELOCITY-X", "Limiter_Velocity_x", "LIMITER", "Limiter value of the x-velocity");
    AddVolumeOutput("LIMITER_VELOCITY-Y", "Limiter_Velocity_y", "LIMITER", "Limiter value of the y-velocity");
    if (nDim == 3) {
      AddVolumeOutput("LIMITER_VELOCITY-Z", "Limiter_Velocity_z", "LIMITER", "Limiter value of the z-velocity");
    }
    AddVolumeOutput("LIMITER_PRESSURE", "Limiter_Pressure", "LIMITER", "Limiter value of the pressure");
    AddVolumeOutput("LIMITER_DENSITY", "Limiter_Density", "LIMITER", "Limiter value of the density");
    AddVolumeOutput("LIMITER_ENTHALPY", "Limiter_Enthalpy", "LIMITER", "Limiter value of the enthalpy");
  }

  SetVolumeOutputFieldsScalarLimiter(config);

  SetVolumeOutputFieldsScalarSource(config);

  SetVolumeOutputFieldsScalarLookup(config);

  SetVolumeOutputFieldsScalarMisc(config);

  if (config->GetCompute_Metric()) {
    AddVolumeOutput("METRIC_XX", "Metric_xx", "SOLUTION", "xx component of mesh adaptation metric tensor");
    AddVolumeOutput("METRIC_XY", "Metric_xy", "SOLUTION", "xy component of mesh adaptation metric tensor");
    AddVolumeOutput("METRIC_YY", "Metric_yy", "SOLUTION", "yy component of mesh adaptation metric tensor");
    if (nDim == 3) {
      AddVolumeOutput("METRIC_XZ", "Metric_xz", "SOLUTION", "xz component of mesh adaptation metric tensor");
      AddVolumeOutput("METRIC_YZ", "Metric_yz", "SOLUTION", "yz component of mesh adaptation metric tensor");
      AddVolumeOutput("METRIC_ZZ", "Metric_zz", "SOLUTION", "zz component of mesh adaptation metric tensor");
    }
  }

  // Roe Low Dissipation
  if (config->GetKind_RoeLowDiss() != NO_ROELOWDISS) {
    AddVolumeOutput("ROE_DISSIPATION", "Roe_Dissipation", "ROE_DISSIPATION", "Value of the Roe dissipation");
  }

  AddCommonFVMOutputs(config);

  if (config->GetTime_Domain()) {
    SetTimeAveragedFields();
  }
}

void CFlowCompOutput::LoadVolumeData(CConfig *config, CGeometry *geometry, CSolver **solver, unsigned long iPoint){
  const auto* Node_Flow = solver[FLOW_SOL]->GetNodes();
  const bool detonation = config->GetKind_Solver() == MAIN_SOLVER::DETONATION_EULER;
  const auto* Detonation_Node_Flow = detonation ? static_cast<const CDetonationEulerVariable*>(Node_Flow) : nullptr;
  const auto* Detonation_Solver = detonation ? static_cast<const CDetonationEulerSolver*>(solver[FLOW_SOL]) : nullptr;
  auto* Node_Geo  = geometry->nodes;

  LoadCoordinates(Node_Geo->GetCoord(iPoint), iPoint);

  if (detonation) {
    const auto nSpecies = Detonation_Node_Flow->GetNumSpecies();
    SetVolumeOutputValue("DENSITY", iPoint, Detonation_Node_Flow->GetDensity(iPoint));
    SetVolumeOutputValue("MOMENTUM-X", iPoint, Node_Flow->GetSolution(iPoint, nSpecies + 0));
    SetVolumeOutputValue("MOMENTUM-Y", iPoint, Node_Flow->GetSolution(iPoint, nSpecies + 1));
    if (nDim == 3) {
      SetVolumeOutputValue("MOMENTUM-Z", iPoint, Node_Flow->GetSolution(iPoint, nSpecies + 2));
      SetVolumeOutputValue("ENERGY", iPoint, Node_Flow->GetSolution(iPoint, nSpecies + 3));
    } else {
      SetVolumeOutputValue("ENERGY", iPoint, Node_Flow->GetSolution(iPoint, nSpecies + 2));
    }
  } else {
    SetVolumeOutputValue("DENSITY",    iPoint, Node_Flow->GetSolution(iPoint, 0));
    SetVolumeOutputValue("MOMENTUM-X", iPoint, Node_Flow->GetSolution(iPoint, 1));
    SetVolumeOutputValue("MOMENTUM-Y", iPoint, Node_Flow->GetSolution(iPoint, 2));
    if (nDim == 3){
      SetVolumeOutputValue("MOMENTUM-Z", iPoint, Node_Flow->GetSolution(iPoint, 3));
      SetVolumeOutputValue("ENERGY",     iPoint, Node_Flow->GetSolution(iPoint, 4));
    } else {
      SetVolumeOutputValue("ENERGY",     iPoint, Node_Flow->GetSolution(iPoint, 3));
    }
  }

  if (gridMovement){
    SetVolumeOutputValue("GRID_VELOCITY-X", iPoint, Node_Geo->GetGridVel(iPoint)[0]);
    SetVolumeOutputValue("GRID_VELOCITY-Y", iPoint, Node_Geo->GetGridVel(iPoint)[1]);
    if (nDim == 3)
      SetVolumeOutputValue("GRID_VELOCITY-Z", iPoint, Node_Geo->GetGridVel(iPoint)[2]);
  }

  SetVolumeOutputValue("PRESSURE", iPoint, Node_Flow->GetPressure(iPoint));
  SetVolumeOutputValue("TEMPERATURE", iPoint, Node_Flow->GetTemperature(iPoint));
  SetVolumeOutputValue("MACH", iPoint, sqrt(Node_Flow->GetVelocity2(iPoint))/Node_Flow->GetSoundSpeed(iPoint));
  SetVolumeOutputValue("SOUND_SPEED", iPoint, Node_Flow->GetSoundSpeed(iPoint));

  const su2double factor = solver[FLOW_SOL]->GetReferenceDynamicPressure();
  SetVolumeOutputValue("PRESSURE_COEFF", iPoint, (Node_Flow->GetPressure(iPoint) - solver[FLOW_SOL]->GetPressure_Inf())/factor);
  SetVolumeOutputValue("VELOCITY-X", iPoint, Node_Flow->GetVelocity(iPoint, 0));
  SetVolumeOutputValue("VELOCITY-Y", iPoint, Node_Flow->GetVelocity(iPoint, 1));
  if (nDim == 3){
    SetVolumeOutputValue("VELOCITY-Z", iPoint, Node_Flow->GetVelocity(iPoint, 2));
  }

  if(config->GetKind_FluidModel() == DATADRIVEN_FLUID){
    SetVolumeOutputValue("EXTRAPOLATION", iPoint, Node_Flow->GetDataExtrapolation(iPoint));
    SetVolumeOutputValue("FLUIDMODEL_NEWTONITER", iPoint, Node_Flow->GetNewtonSolverIterations(iPoint));
    SetVolumeOutputValue("ENTROPY", iPoint, Node_Flow->GetEntropy(iPoint));
  }

  if (config->GetKind_Solver() == MAIN_SOLVER::RANS || config->GetKind_Solver() == MAIN_SOLVER::NAVIER_STOKES){
    SetVolumeOutputValue("LAMINAR_VISCOSITY", iPoint, Node_Flow->GetLaminarViscosity(iPoint));
  }

  if (detonation) {
    const auto nSpecies = Detonation_Node_Flow->GetNumSpecies();
    su2double densityResidual = 0.0;
    for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      densityResidual += solver[FLOW_SOL]->LinSysRes(iPoint, iSpecies);
    }
    SetVolumeOutputValue("RES_DENSITY", iPoint, densityResidual);
    SetVolumeOutputValue("RES_MOMENTUM-X", iPoint, solver[FLOW_SOL]->LinSysRes(iPoint, nSpecies + 0));
    SetVolumeOutputValue("RES_MOMENTUM-Y", iPoint, solver[FLOW_SOL]->LinSysRes(iPoint, nSpecies + 1));
    if (nDim == 3) {
      SetVolumeOutputValue("RES_MOMENTUM-Z", iPoint, solver[FLOW_SOL]->LinSysRes(iPoint, nSpecies + 2));
      SetVolumeOutputValue("RES_ENERGY", iPoint, solver[FLOW_SOL]->LinSysRes(iPoint, nSpecies + 3));
    } else {
      SetVolumeOutputValue("RES_ENERGY", iPoint, solver[FLOW_SOL]->LinSysRes(iPoint, nSpecies + 2));
    }
  } else {
    SetVolumeOutputValue("RES_DENSITY", iPoint, solver[FLOW_SOL]->LinSysRes(iPoint, 0));
    SetVolumeOutputValue("RES_MOMENTUM-X", iPoint, solver[FLOW_SOL]->LinSysRes(iPoint, 1));
    SetVolumeOutputValue("RES_MOMENTUM-Y", iPoint, solver[FLOW_SOL]->LinSysRes(iPoint, 2));
    if (nDim == 3){
      SetVolumeOutputValue("RES_MOMENTUM-Z", iPoint, solver[FLOW_SOL]->LinSysRes(iPoint, 3));
      SetVolumeOutputValue("RES_ENERGY", iPoint, solver[FLOW_SOL]->LinSysRes(iPoint, 4));
    } else {
      SetVolumeOutputValue("RES_ENERGY", iPoint, solver[FLOW_SOL]->LinSysRes(iPoint, 3));
    }
  }

  if (config->GetKind_SlopeLimit_Flow() != LIMITER::NONE && config->GetKind_SlopeLimit_Flow() != LIMITER::VAN_ALBADA_EDGE) {
    if (detonation) {
      const auto nSpecies = Detonation_Node_Flow->GetNumSpecies();
      SetVolumeOutputValue("LIMITER_VELOCITY-X", iPoint, Node_Flow->GetLimiter_Primitive(iPoint, nSpecies + 1));
      SetVolumeOutputValue("LIMITER_VELOCITY-Y", iPoint, Node_Flow->GetLimiter_Primitive(iPoint, nSpecies + 2));
      if (nDim == 3){
        SetVolumeOutputValue("LIMITER_VELOCITY-Z", iPoint, Node_Flow->GetLimiter_Primitive(iPoint, nSpecies + 3));
      }
      SetVolumeOutputValue("LIMITER_PRESSURE", iPoint, Node_Flow->GetLimiter_Primitive(iPoint, nSpecies + nDim + 1));
      SetVolumeOutputValue("LIMITER_DENSITY", iPoint, Node_Flow->GetLimiter_Primitive(iPoint, nSpecies + nDim + 2));
      SetVolumeOutputValue("LIMITER_ENTHALPY", iPoint, Node_Flow->GetLimiter_Primitive(iPoint, nSpecies + nDim + 3));
    } else {
      SetVolumeOutputValue("LIMITER_VELOCITY-X", iPoint, Node_Flow->GetLimiter_Primitive(iPoint, 1));
      SetVolumeOutputValue("LIMITER_VELOCITY-Y", iPoint, Node_Flow->GetLimiter_Primitive(iPoint, 2));
      if (nDim == 3){
        SetVolumeOutputValue("LIMITER_VELOCITY-Z", iPoint, Node_Flow->GetLimiter_Primitive(iPoint, 3));
      }
      SetVolumeOutputValue("LIMITER_PRESSURE", iPoint, Node_Flow->GetLimiter_Primitive(iPoint, nDim+1));
      SetVolumeOutputValue("LIMITER_DENSITY", iPoint, Node_Flow->GetLimiter_Primitive(iPoint, nDim+2));
      SetVolumeOutputValue("LIMITER_ENTHALPY", iPoint, Node_Flow->GetLimiter_Primitive(iPoint, nDim+3));
    }
  }

  if (detonation) {
    const auto nSpecies = Detonation_Node_Flow->GetNumSpecies();
    for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      SetVolumeOutputValue("RHOS_" + std::to_string(iSpecies), iPoint, Node_Flow->GetSolution(iPoint, iSpecies));
      SetVolumeOutputValue("MASSFRAC_" + std::to_string(iSpecies), iPoint,
                           Detonation_Node_Flow->GetMassFraction(iPoint, iSpecies));
    }
    SetVolumeOutputValue("QDOT", iPoint, Detonation_Solver->GetChemQdot(iPoint));
    SetVolumeOutputValue("CHEM_FAIL_FLAG", iPoint, Detonation_Solver->GetChemFailFlag(iPoint));
    SetVolumeOutputValue("CHEM_FAIL_COUNT", iPoint, Detonation_Solver->GetChemFailCount(iPoint));
    SetVolumeOutputValue("CHEM_SUBSTEPS", iPoint, Detonation_Solver->GetChemSubstepsLast(iPoint));
  }

  if (config->GetKind_RoeLowDiss() != NO_ROELOWDISS){
    SetVolumeOutputValue("ROE_DISSIPATION", iPoint, Node_Flow->GetRoe_Dissipation(iPoint));
  }

  if (config->GetCompute_Metric()) {
    SetVolumeOutputValue("METRIC_XX", iPoint, Node_Flow->GetMetric(iPoint, 0));
    SetVolumeOutputValue("METRIC_XY", iPoint, Node_Flow->GetMetric(iPoint, 1));
    SetVolumeOutputValue("METRIC_YY", iPoint, Node_Flow->GetMetric(iPoint, 2));
    if (nDim == 3) {
      SetVolumeOutputValue("METRIC_XZ", iPoint, Node_Flow->GetMetric(iPoint, 3));
      SetVolumeOutputValue("METRIC_YZ", iPoint, Node_Flow->GetMetric(iPoint, 4));
      SetVolumeOutputValue("METRIC_ZZ", iPoint, Node_Flow->GetMetric(iPoint, 5));
    }
  }

  LoadVolumeDataScalar(config, solver, geometry, iPoint);

  LoadCommonFVMOutputs(config, geometry, iPoint);

  if (config->GetTime_Domain()) {
    LoadTimeAveragedData(iPoint, Node_Flow);
  }
}

void CFlowCompOutput::LoadHistoryData(CConfig *config, CGeometry *geometry, CSolver **solver)  {

  CSolver* flow_solver = solver[FLOW_SOL];
  CSolver* mesh_solver = solver[MESH_SOL];
  const bool detonation = config->GetKind_Solver() == MAIN_SOLVER::DETONATION_EULER;
  const auto* detonation_solver = detonation ? static_cast<const CDetonationEulerSolver*>(flow_solver) : nullptr;
  if (detonation) {
    const auto nSpecies = detonation_solver->GetNumSpecies();
    su2double densityRMS = 0.0, densityMax = 0.0;
    for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      densityRMS += pow(flow_solver->GetRes_RMS(iSpecies), 2);
      densityMax = max(densityMax, flow_solver->GetRes_Max(iSpecies));
    }
    densityRMS = sqrt(densityRMS);

    SetHistoryOutputValue("RMS_DENSITY", log10(densityRMS));
    SetHistoryOutputValue("RMS_MOMENTUM-X", log10(flow_solver->GetRes_RMS(nSpecies + 0)));
    SetHistoryOutputValue("RMS_MOMENTUM-Y", log10(flow_solver->GetRes_RMS(nSpecies + 1)));
    if (nDim == 2) {
      SetHistoryOutputValue("RMS_ENERGY", log10(flow_solver->GetRes_RMS(nSpecies + 2)));
    } else {
      SetHistoryOutputValue("RMS_MOMENTUM-Z", log10(flow_solver->GetRes_RMS(nSpecies + 2)));
      SetHistoryOutputValue("RMS_ENERGY", log10(flow_solver->GetRes_RMS(nSpecies + 3)));
    }

    SetHistoryOutputValue("MAX_DENSITY", log10(densityMax));
    SetHistoryOutputValue("MAX_MOMENTUM-X", log10(flow_solver->GetRes_Max(nSpecies + 0)));
    SetHistoryOutputValue("MAX_MOMENTUM-Y", log10(flow_solver->GetRes_Max(nSpecies + 1)));
    if (nDim == 2) {
      SetHistoryOutputValue("MAX_ENERGY", log10(flow_solver->GetRes_Max(nSpecies + 2)));
    } else {
      SetHistoryOutputValue("MAX_MOMENTUM-Z", log10(flow_solver->GetRes_Max(nSpecies + 2)));
      SetHistoryOutputValue("MAX_ENERGY", log10(flow_solver->GetRes_Max(nSpecies + 3)));
    }

    if (multiZone) {
      su2double densityBGS = 0.0;
      for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
        densityBGS += pow(flow_solver->GetRes_BGS(iSpecies), 2);
      }
      densityBGS = sqrt(densityBGS);

      SetHistoryOutputValue("BGS_DENSITY", log10(densityBGS));
      SetHistoryOutputValue("BGS_MOMENTUM-X", log10(flow_solver->GetRes_BGS(nSpecies + 0)));
      SetHistoryOutputValue("BGS_MOMENTUM-Y", log10(flow_solver->GetRes_BGS(nSpecies + 1)));
      if (nDim == 2) {
        SetHistoryOutputValue("BGS_ENERGY", log10(flow_solver->GetRes_BGS(nSpecies + 2)));
      } else {
        SetHistoryOutputValue("BGS_MOMENTUM-Z", log10(flow_solver->GetRes_BGS(nSpecies + 2)));
        SetHistoryOutputValue("BGS_ENERGY", log10(flow_solver->GetRes_BGS(nSpecies + 3)));
      }
    }
  } else {
    SetHistoryOutputValue("RMS_DENSITY", log10(flow_solver->GetRes_RMS(0)));
    SetHistoryOutputValue("RMS_MOMENTUM-X", log10(flow_solver->GetRes_RMS(1)));
    SetHistoryOutputValue("RMS_MOMENTUM-Y", log10(flow_solver->GetRes_RMS(2)));
    if (nDim == 2)
      SetHistoryOutputValue("RMS_ENERGY", log10(flow_solver->GetRes_RMS(3)));
    else {
      SetHistoryOutputValue("RMS_MOMENTUM-Z", log10(flow_solver->GetRes_RMS(3)));
      SetHistoryOutputValue("RMS_ENERGY", log10(flow_solver->GetRes_RMS(4)));
    }
    SetHistoryOutputValue("MAX_DENSITY", log10(flow_solver->GetRes_Max(0)));
    SetHistoryOutputValue("MAX_MOMENTUM-X", log10(flow_solver->GetRes_Max(1)));
    SetHistoryOutputValue("MAX_MOMENTUM-Y", log10(flow_solver->GetRes_Max(2)));
    if (nDim == 2)
      SetHistoryOutputValue("MAX_ENERGY", log10(flow_solver->GetRes_Max(3)));
    else {
      SetHistoryOutputValue("MAX_MOMENTUM-Z", log10(flow_solver->GetRes_Max(3)));
      SetHistoryOutputValue("MAX_ENERGY", log10(flow_solver->GetRes_Max(4)));
    }
    if (multiZone){
      SetHistoryOutputValue("BGS_DENSITY", log10(flow_solver->GetRes_BGS(0)));
      SetHistoryOutputValue("BGS_MOMENTUM-X", log10(flow_solver->GetRes_BGS(1)));
      SetHistoryOutputValue("BGS_MOMENTUM-Y", log10(flow_solver->GetRes_BGS(2)));
      if (nDim == 2)
        SetHistoryOutputValue("BGS_ENERGY", log10(flow_solver->GetRes_BGS(3)));
      else {
        SetHistoryOutputValue("BGS_MOMENTUM-Z", log10(flow_solver->GetRes_BGS(3)));
        SetHistoryOutputValue("BGS_ENERGY", log10(flow_solver->GetRes_BGS(4)));
      }
    }
  }

  SetHistoryOutputValue("MIN_DELTA_TIME", flow_solver->GetMin_Delta_Time());
  SetHistoryOutputValue("MAX_DELTA_TIME", flow_solver->GetMax_Delta_Time());

  SetHistoryOutputValue("MIN_CFL", flow_solver->GetMin_CFL_Local());
  SetHistoryOutputValue("MAX_CFL", flow_solver->GetMax_CFL_Local());
  SetHistoryOutputValue("AVG_CFL", flow_solver->GetAvg_CFL_Local());

  if (detonation) {
    const auto* detonation_nodes = static_cast<const CDetonationEulerVariable*>(flow_solver->GetNodes());
    const auto nSpecies = detonation_solver->GetNumSpecies();

    su2double localSumTemperature = 0.0;
    su2double localSumPressure = 0.0;
    su2double localSumDensity = 0.0;
    su2double localSumInternalEnergy = 0.0;
    su2double localSumQdot = 0.0;
    su2double localMinTemperature = std::numeric_limits<su2double>::max();
    su2double localMaxTemperature = std::numeric_limits<su2double>::lowest();
    su2double localMinPressure = std::numeric_limits<su2double>::max();
    su2double localMaxPressure = std::numeric_limits<su2double>::lowest();
    su2double localMinDensity = std::numeric_limits<su2double>::max();
    su2double localMaxDensity = std::numeric_limits<su2double>::lowest();
    std::vector<su2double> localMinMomentum(nDim, std::numeric_limits<su2double>::max());
    std::vector<su2double> localMaxMomentum(nDim, std::numeric_limits<su2double>::lowest());
    su2double localMinEnergy = std::numeric_limits<su2double>::max();
    su2double localMaxEnergy = std::numeric_limits<su2double>::lowest();
    su2double localMinMassFraction = std::numeric_limits<su2double>::max();
    su2double localMaxSumYError = 0.0;
    std::vector<su2double> localSumY(nSpecies, 0.0);
    std::vector<su2double> localMinY(nSpecies, std::numeric_limits<su2double>::max());
    std::vector<su2double> localMaxY(nSpecies, std::numeric_limits<su2double>::lowest());
    std::vector<su2double> localMinRhoS(nSpecies, std::numeric_limits<su2double>::max());
    std::vector<su2double> localMaxRhoS(nSpecies, std::numeric_limits<su2double>::lowest());
    unsigned long localDomainPoints = 0;

    for (unsigned long iPoint = 0; iPoint < geometry->GetnPoint(); ++iPoint) {
      if (!geometry->nodes->GetDomain(iPoint)) continue;

      localDomainPoints++;

      const su2double temperature = detonation_nodes->GetTemperature(iPoint);
      const su2double pressure = detonation_nodes->GetPressure(iPoint);
      const su2double density = detonation_nodes->GetDensity(iPoint);
      const su2double internalEnergy = detonation_nodes->GetEnergy(iPoint) - 0.5 * detonation_nodes->GetVelocity2(iPoint);

      localSumTemperature += temperature;
      localSumPressure += pressure;
      localSumDensity += density;
      localSumInternalEnergy += internalEnergy;
      localSumQdot += detonation_solver->GetChemQdot(iPoint);

      localMinTemperature = min(localMinTemperature, temperature);
      localMaxTemperature = max(localMaxTemperature, temperature);
      localMinPressure = min(localMinPressure, pressure);
      localMaxPressure = max(localMaxPressure, pressure);
      localMinDensity = min(localMinDensity, density);
      localMaxDensity = max(localMaxDensity, density);
      for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
        const su2double momentum = flow_solver->GetNodes()->GetSolution(iPoint, nSpecies + iDim);
        localMinMomentum[iDim] = min(localMinMomentum[iDim], momentum);
        localMaxMomentum[iDim] = max(localMaxMomentum[iDim], momentum);
      }
      const su2double totalEnergy = flow_solver->GetNodes()->GetSolution(iPoint, nSpecies + nDim);
      localMinEnergy = min(localMinEnergy, totalEnergy);
      localMaxEnergy = max(localMaxEnergy, totalEnergy);

      su2double sumY = 0.0;
      for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
        const su2double y = detonation_nodes->GetMassFraction(iPoint, iSpecies);
        const su2double rhoS = flow_solver->GetNodes()->GetSolution(iPoint, iSpecies);
        localSumY[iSpecies] += y;
        localMinY[iSpecies] = min(localMinY[iSpecies], y);
        localMaxY[iSpecies] = max(localMaxY[iSpecies], y);
        localMinRhoS[iSpecies] = min(localMinRhoS[iSpecies], rhoS);
        localMaxRhoS[iSpecies] = max(localMaxRhoS[iSpecies], rhoS);
        localMinMassFraction = min(localMinMassFraction, y);
        sumY += y;
      }
      localMaxSumYError = max(localMaxSumYError, fabs(sumY - 1.0));
    }

    unsigned long globalDomainPoints = localDomainPoints;
    su2double globalSumTemperature = localSumTemperature;
    su2double globalSumPressure = localSumPressure;
    su2double globalSumDensity = localSumDensity;
    su2double globalSumInternalEnergy = localSumInternalEnergy;
    su2double globalSumQdot = localSumQdot;
    su2double globalMinTemperature = localMinTemperature;
    su2double globalMaxTemperature = localMaxTemperature;
    su2double globalMinPressure = localMinPressure;
    su2double globalMaxPressure = localMaxPressure;
    su2double globalMinDensity = localMinDensity;
    su2double globalMaxDensity = localMaxDensity;
    std::vector<su2double> globalMinMomentum = localMinMomentum;
    std::vector<su2double> globalMaxMomentum = localMaxMomentum;
    su2double globalMinEnergy = localMinEnergy;
    su2double globalMaxEnergy = localMaxEnergy;
    su2double globalMinMassFraction = localMinMassFraction;
    su2double globalMaxSumYError = localMaxSumYError;
    std::vector<su2double> globalSumY = localSumY;
    std::vector<su2double> globalMinY = localMinY;
    std::vector<su2double> globalMaxY = localMaxY;
    std::vector<su2double> globalMinRhoS = localMinRhoS;
    std::vector<su2double> globalMaxRhoS = localMaxRhoS;

    if (config->GetComm_Level() == COMM_FULL) {
      SU2_MPI::Allreduce(&localDomainPoints, &globalDomainPoints, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
      SU2_MPI::Allreduce(&localSumTemperature, &globalSumTemperature, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
      SU2_MPI::Allreduce(&localSumPressure, &globalSumPressure, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
      SU2_MPI::Allreduce(&localSumDensity, &globalSumDensity, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
      SU2_MPI::Allreduce(&localSumInternalEnergy, &globalSumInternalEnergy, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
      SU2_MPI::Allreduce(&localSumQdot, &globalSumQdot, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
      SU2_MPI::Allreduce(&localMinTemperature, &globalMinTemperature, 1, MPI_DOUBLE, MPI_MIN, SU2_MPI::GetComm());
      SU2_MPI::Allreduce(&localMaxTemperature, &globalMaxTemperature, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());
      SU2_MPI::Allreduce(&localMinPressure, &globalMinPressure, 1, MPI_DOUBLE, MPI_MIN, SU2_MPI::GetComm());
      SU2_MPI::Allreduce(&localMaxPressure, &globalMaxPressure, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());
      SU2_MPI::Allreduce(&localMinDensity, &globalMinDensity, 1, MPI_DOUBLE, MPI_MIN, SU2_MPI::GetComm());
      SU2_MPI::Allreduce(&localMaxDensity, &globalMaxDensity, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());
      SU2_MPI::Allreduce(localMinMomentum.data(), globalMinMomentum.data(), nDim, MPI_DOUBLE, MPI_MIN, SU2_MPI::GetComm());
      SU2_MPI::Allreduce(localMaxMomentum.data(), globalMaxMomentum.data(), nDim, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());
      SU2_MPI::Allreduce(&localMinEnergy, &globalMinEnergy, 1, MPI_DOUBLE, MPI_MIN, SU2_MPI::GetComm());
      SU2_MPI::Allreduce(&localMaxEnergy, &globalMaxEnergy, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());
      SU2_MPI::Allreduce(&localMinMassFraction, &globalMinMassFraction, 1, MPI_DOUBLE, MPI_MIN, SU2_MPI::GetComm());
      SU2_MPI::Allreduce(&localMaxSumYError, &globalMaxSumYError, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());
      SU2_MPI::Allreduce(localSumY.data(), globalSumY.data(), nSpecies, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
      SU2_MPI::Allreduce(localMinY.data(), globalMinY.data(), nSpecies, MPI_DOUBLE, MPI_MIN, SU2_MPI::GetComm());
      SU2_MPI::Allreduce(localMaxY.data(), globalMaxY.data(), nSpecies, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());
      SU2_MPI::Allreduce(localMinRhoS.data(), globalMinRhoS.data(), nSpecies, MPI_DOUBLE, MPI_MIN, SU2_MPI::GetComm());
      SU2_MPI::Allreduce(localMaxRhoS.data(), globalMaxRhoS.data(), nSpecies, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());
    }

    const su2double invDomainPoints = (globalDomainPoints > 0) ? 1.0 / static_cast<su2double>(globalDomainPoints) : 0.0;
    unsigned long localHLLCFallbackCount = detonation_solver->GetHLLCFallbackCount();
    unsigned long localHLLCNegativeStateCount = detonation_solver->GetHLLCNegativeStateCount();
    unsigned long globalHLLCFallbackCount = localHLLCFallbackCount;
    unsigned long globalHLLCNegativeStateCount = localHLLCNegativeStateCount;
    su2double localHLLCMinPressure = detonation_solver->GetHLLCMinPressure();
    su2double localHLLCMinDensity = detonation_solver->GetHLLCMinDensity();
    su2double globalHLLCMinPressure = localHLLCMinPressure;
    su2double globalHLLCMinDensity = localHLLCMinDensity;
    if (config->GetComm_Level() == COMM_FULL) {
      SU2_MPI::Allreduce(&localHLLCFallbackCount, &globalHLLCFallbackCount, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
      SU2_MPI::Allreduce(&localHLLCNegativeStateCount, &globalHLLCNegativeStateCount, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
      SU2_MPI::Allreduce(&localHLLCMinPressure, &globalHLLCMinPressure, 1, MPI_DOUBLE, MPI_MIN, SU2_MPI::GetComm());
      SU2_MPI::Allreduce(&localHLLCMinDensity, &globalHLLCMinDensity, 1, MPI_DOUBLE, MPI_MIN, SU2_MPI::GetComm());
    }
    SetHistoryOutputValue("CHEM_FAIL_CELLS", detonation_solver->GetChemFailCells());
    SetHistoryOutputValue("CHEM_SUBSTEPS_MAX", detonation_solver->GetChemSubstepsMax());
    SetHistoryOutputValue("CHEM_E_ABS_MAX", detonation_solver->GetChemEnergyAbsMax());
    SetHistoryOutputValue("CHEM_E_REL_MAX", detonation_solver->GetChemEnergyRelMax());
    SetHistoryOutputValue("HLLC_FALLBACK_COUNT", globalHLLCFallbackCount);
    SetHistoryOutputValue("HLLC_NEGATIVE_STATE_COUNT", globalHLLCNegativeStateCount);
    SetHistoryOutputValue("HLLC_MIN_PRESSURE", globalHLLCMinPressure);
    SetHistoryOutputValue("HLLC_MIN_DENSITY", globalHLLCMinDensity);
    SetHistoryOutputValue("DET_AVG_TEMPERATURE", globalSumTemperature * invDomainPoints);
    SetHistoryOutputValue("DET_AVG_PRESSURE", globalSumPressure * invDomainPoints);
    SetHistoryOutputValue("DET_AVG_DENSITY", globalSumDensity * invDomainPoints);
    SetHistoryOutputValue("DET_AVG_INTERNAL_ENERGY", globalSumInternalEnergy * invDomainPoints);
    SetHistoryOutputValue("DET_AVG_QDOT", globalSumQdot * invDomainPoints);
    SetHistoryOutputValue("DET_SPAN_TEMPERATURE", globalMaxTemperature - globalMinTemperature);
    SetHistoryOutputValue("DET_SPAN_PRESSURE", globalMaxPressure - globalMinPressure);
    SetHistoryOutputValue("DET_SPAN_DENSITY", globalMaxDensity - globalMinDensity);
    SetHistoryOutputValue("DET_SPAN_MOMENTUM-X", globalMaxMomentum[0] - globalMinMomentum[0]);
    SetHistoryOutputValue("DET_SPAN_MOMENTUM-Y", globalMaxMomentum[1] - globalMinMomentum[1]);
    if (nDim == 3) {
      SetHistoryOutputValue("DET_SPAN_MOMENTUM-Z", globalMaxMomentum[2] - globalMinMomentum[2]);
    }
    SetHistoryOutputValue("DET_SPAN_ENERGY", globalMaxEnergy - globalMinEnergy);
    SetHistoryOutputValue("DET_MIN_MASS_FRACTION", globalMinMassFraction);
    SetHistoryOutputValue("DET_MAX_SUM_Y_ERROR", globalMaxSumYError);
    for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      const auto suffix = std::to_string(iSpecies);
      SetHistoryOutputValue("DET_AVG_Y_" + suffix, globalSumY[iSpecies] * invDomainPoints);
      SetHistoryOutputValue("DET_SPAN_Y_" + suffix, globalMaxY[iSpecies] - globalMinY[iSpecies]);
      SetHistoryOutputValue("DET_SPAN_RHOS_" + suffix, globalMaxRhoS[iSpecies] - globalMinRhoS[iSpecies]);
    }
  }

  SetHistoryOutputValue("LINSOL_ITER", flow_solver->GetIterLinSolver());
  SetHistoryOutputValue("LINSOL_RESIDUAL", log10(flow_solver->GetResLinSolver()));

  if (config->GetDeform_Mesh()){
    SetHistoryOutputValue("DEFORM_MIN_VOLUME", mesh_solver->GetMinimum_Volume());
    SetHistoryOutputValue("DEFORM_MAX_VOLUME", mesh_solver->GetMaximum_Volume());
    SetHistoryOutputValue("DEFORM_ITER", mesh_solver->GetIterLinSolver());
    SetHistoryOutputValue("DEFORM_RESIDUAL", log10(mesh_solver->GetResLinSolver()));
  }

  if(config->GetFixed_CL_Mode()){
    SetHistoryOutputValue("DELTA_CL", fabs(flow_solver->GetTotal_CL() - config->GetTarget_CL()));
    SetHistoryOutputValue("PREV_AOA", flow_solver->GetPrevious_AoA());
    SetHistoryOutputValue("CHANGE_IN_AOA", config->GetAoA()-flow_solver->GetPrevious_AoA());
    SetHistoryOutputValue("CL_DRIVER_COMMAND", flow_solver->GetAoA_inc());
  }

  LoadHistoryDataScalar(config, solver);

  /*--- Set the analyse surface history values --- */

  SetAnalyzeSurface(solver, geometry, config, false);

  /*--- Set aerodynamic coefficients --- */

  SetAerodynamicCoefficients(config, flow_solver);

  if (config->GetViscous()) {
    SetHistoryOutputValue("BUFFET", flow_solver->GetTotal_Buffet_Metric());
  }

  SetHeatCoefficients(config, flow_solver);

  /*--- Set rotating frame coefficients --- */

  SetRotatingFrameCoefficients(flow_solver);

  /*--- Set Cp diff fields ---*/

  SetCpInverseDesign(flow_solver, geometry, config);

  /*--- Set nearfield diff fields ---*/

  if (config->GetEquivArea()) SetNearfieldInverseDesign(flow_solver, geometry, config);

  /*--- Keep this as last, since it uses the history values that were set. ---*/

  SetCustomOutputs(solver, geometry, config);

  SetCustomAndComboObjectives(FLOW_SOL, config, solver);
}

bool CFlowCompOutput::SetInitResiduals(const CConfig *config){

  return (config->GetTime_Marching() != TIME_MARCHING::STEADY && (curInnerIter == 0))||
         (config->GetTime_Marching() == TIME_MARCHING::STEADY && (curInnerIter < 2));

}

void CFlowCompOutput::SetAdditionalScreenOutput(const CConfig *config){

  if (config->GetFixed_CL_Mode()){
    SetFixedCLScreenOutput(config);
  }
}

bool CFlowCompOutput::WriteHistoryFileOutput(const CConfig *config) {
  return !config->GetFinite_Difference_Mode() && COutput::WriteHistoryFileOutput(config);
}

void CFlowCompOutput::SetTurboPerformance_Output(std::shared_ptr<CTurboOutput> TurboPerf,
                                  CConfig *config,
                                  unsigned long TimeIter,
                                  unsigned long OuterIter,
                                  unsigned long InnerIter) {

  curTimeIter  = TimeIter;
  curAbsTimeIter = TimeIter - config->GetRestart_Iter();
  curOuterIter = OuterIter;
  curInnerIter = InnerIter;
  stringstream TurboInOutTable, TurboPerfTable;

  auto BladePerformance = TurboPerf->GetBladesPerformances();

  /*-- Table for Turbomachinery Performance Values --*/
  PrintingToolbox::CTablePrinter TurboInOut(&TurboInOutTable);

  TurboInOutTable<<"-- Turbomachinery inlet and outlet property Summary:"<<endl;
  TurboInOut.AddColumn("Properties", 25);
  TurboInOut.AddColumn("Inlet", 25);
  TurboInOut.AddColumn("Outlet", 25);
  TurboInOut.SetAlign(PrintingToolbox::CTablePrinter::RIGHT);
  TurboInOut.PrintHeader();

  for (unsigned short iZone = 0; iZone <= config->GetnZone()-1; iZone++) {
    auto nSpan = config->GetnSpan_iZones(iZone);
    const auto& BladePerf = BladePerformance.at(iZone).at(nSpan);

    TurboInOut<<" BLADE ROW INDEX "<<iZone <<"";
    TurboInOut.PrintFooter();
    // TODO: Blade Wise Printing
    TurboInOut << "Entropy " << BladePerf->GetInletState().GetEntropy() << BladePerf->GetOutletState().GetEntropy();
    TurboInOut << "Total Enthalpy " << BladePerf->GetInletState().GetTotalEnthalpy() << BladePerf->GetOutletState().GetTotalEnthalpy();
    TurboInOut << "Total Pressure " << BladePerf->GetInletState().GetTotalPressure() << BladePerf->GetOutletState().GetTotalPressure();
    TurboInOut << "Pressure " << BladePerf->GetInletState().GetPressure() << BladePerf->GetOutletState().GetPressure();
    TurboInOut << "Density " << BladePerf->GetInletState().GetDensity() << BladePerf->GetOutletState().GetDensity();
    TurboInOut << "Normal Velocity " << BladePerf->GetInletState().GetVelocity()[0] << BladePerf->GetOutletState().GetVelocity()[0];
    TurboInOut << "Tangential Velocity " << BladePerf->GetInletState().GetVelocity()[1] << BladePerf->GetOutletState().GetVelocity()[1];
    TurboInOut << "Mass Flow " << BladePerf->GetInletState().GetMassFlow() << BladePerf->GetOutletState().GetMassFlow();
    TurboInOut << "Mach " << BladePerf->GetInletState().GetMachValue() << BladePerf->GetOutletState().GetMachValue();
    TurboInOut << "Abs Flow Angle " << BladePerf->GetInletState().GetAbsFlowAngle()*180/PI_NUMBER << BladePerf->GetOutletState().GetAbsFlowAngle()*180/PI_NUMBER;
    TurboInOut.PrintFooter();
  }
  cout<<TurboInOutTable.str();
}

void CFlowCompOutput::SetTurboMultiZonePerformance_Output(std::shared_ptr<CTurbomachineryStagePerformance> TurboStagePerf, std::shared_ptr<CTurboOutput> TurboPerf, CConfig *config) {

  stringstream TurboMZPerf;

  PrintingToolbox::CTablePrinter TurboInOut(&TurboMZPerf);

  /*--- Print header for the stage performance computation ---*/
  TurboMZPerf<<"-- Turbomachinery Stage Performance --"<<endl;
  TurboInOut.AddColumn("Index", 13);
  TurboInOut.AddColumn(" Sgen    (%)", 13);
  TurboInOut.AddColumn(" Work (J/kg)", 13);
  TurboInOut.AddColumn(" Efi ts  (%)", 13);
  TurboInOut.AddColumn(" Efi tt  (%)", 13);
  TurboInOut.AddColumn(" PR ts   (-)", 13);
  TurboInOut.AddColumn(" PR tt   (-)", 13);
  TurboInOut.SetAlign(PrintingToolbox::CTablePrinter::RIGHT);
  TurboInOut.PrintHeader();

  /*--- Print Machine Performance (In future also add if the performance is TURBINE or COMPRESSOR) ---*/
  TurboInOut<<"MACHINE"<<TurboStagePerf->GetNormEntropyGen()*100
                        <<TurboStagePerf->GetEulerianWork()
                        <<TurboStagePerf->GetTotalStaticEfficiency()*100
                        <<TurboStagePerf->GetTotalTotalEfficiency()*100
                        <<TurboStagePerf->GetTotalStaticPressureRatio()
                        <<TurboStagePerf->GetTotalTotalPressureRatio();
  TurboInOut.PrintFooter();
  cout<<TurboMZPerf.str();

}

void CFlowCompOutput::LoadTurboHistoryData(std::shared_ptr<CTurbomachineryStagePerformance> TurboStagePerf, std::shared_ptr<CTurboOutput> TurboPerf, CConfig *config) {
  auto BladePerformance = TurboPerf->GetBladesPerformances();
  for (unsigned short iZone = 0; iZone <= config->GetnZone()-1; iZone++) {
    auto nSpan = config->GetnSpan_iZones(iZone);
    const auto& BladePerf = BladePerformance.at(iZone).at(nSpan);

    stringstream tag;
    tag << iZone + 1;

    SetHistoryOutputValue("EntropyIn_" + tag.str(), BladePerf->GetInletState().GetEntropy());
    SetHistoryOutputValue("EntropyOut_" + tag.str(), BladePerf->GetOutletState().GetEntropy());
    SetHistoryOutputValue("TotalEntahalpyIn_" + tag.str(), BladePerf->GetInletState().GetTotalEnthalpy());
    SetHistoryOutputValue("TotalEnthalpyOut_" + tag.str(), BladePerf->GetOutletState().GetTotalEnthalpy());
    SetHistoryOutputValue("TotalPressureIn_" + tag.str(), BladePerf->GetInletState().GetTotalPressure());
    SetHistoryOutputValue("TotalPressureOut_" + tag.str(), BladePerf->GetOutletState().GetTotalPressure());
    SetHistoryOutputValue("PressureIn_" + tag.str(), BladePerf->GetInletState().GetPressure());
    SetHistoryOutputValue("PressureOut_" + tag.str(), BladePerf->GetOutletState().GetPressure());
    SetHistoryOutputValue("TotalTemperatureIn_" + tag.str(), BladePerf->GetInletState().GetTotalTemperature());
    SetHistoryOutputValue("TotalTemperatureOut_" + tag.str(), BladePerf->GetOutletState().GetTotalTemperature());
    SetHistoryOutputValue("TemperatureIn_" + tag.str(), BladePerf->GetInletState().GetTemperature());
    SetHistoryOutputValue("TemperatureOut_" + tag.str(), BladePerf->GetOutletState().GetTemperature());
    SetHistoryOutputValue("DensityIn_" + tag.str(), BladePerf->GetInletState().GetDensity());
    SetHistoryOutputValue("DensityOut_" + tag.str(), BladePerf->GetOutletState().GetDensity());
    SetHistoryOutputValue("NormalVelocityIn_" + tag.str(), BladePerf->GetInletState().GetVelocity()[0]);
    SetHistoryOutputValue("NormalVelocityOut_" + tag.str(), BladePerf->GetOutletState().GetVelocity()[0]);
    SetHistoryOutputValue("TangentialVelocityIn_" + tag.str(), BladePerf->GetInletState().GetVelocity()[1]);
    SetHistoryOutputValue("TangentialVelocityOut_" + tag.str(), BladePerf->GetOutletState().GetVelocity()[1]);
    SetHistoryOutputValue("MassFlowIn_" + tag.str(), BladePerf->GetInletState().GetMassFlow());
    SetHistoryOutputValue("MassFlowOut_" + tag.str(), BladePerf->GetOutletState().GetMassFlow());
    SetHistoryOutputValue("MachIn_" + tag.str(), BladePerf->GetInletState().GetMachValue());
    SetHistoryOutputValue("MachOut_" + tag.str(), BladePerf->GetOutletState().GetMachValue());
    SetHistoryOutputValue("AbsFlowAngleIn_" + tag.str(), BladePerf->GetInletState().GetAbsFlowAngle()*180/PI_NUMBER);
    SetHistoryOutputValue("AbsFlowAngleOut_" + tag.str(), BladePerf->GetOutletState().GetAbsFlowAngle()*180/PI_NUMBER);
    SetHistoryOutputValue("KineticEnergyLoss_" + tag.str(), BladePerf->GetKineticEnergyLoss());
    SetHistoryOutputValue("TotPressureLoss_" + tag.str(), BladePerf->GetTotalPressureLoss());
  }
  SetHistoryOutputValue("EntropyGeneration", TurboStagePerf->GetNormEntropyGen()*100);
  SetHistoryOutputValue("EulerianWork", TurboStagePerf->GetEulerianWork());
  SetHistoryOutputValue("TotalStaticEfficiency", TurboStagePerf->GetTotalStaticEfficiency()*100);
  SetHistoryOutputValue("TotalTotalEfficiency", TurboStagePerf->GetTotalTotalEfficiency()*100);
  SetHistoryOutputValue("PressureRatioTS", TurboStagePerf->GetTotalStaticPressureRatio());
  SetHistoryOutputValue("PressureRatioTT", TurboStagePerf->GetTotalTotalPressureRatio());
  SetHistoryOutputValue("KineticEnergyLoss_Stage", TurboStagePerf->GetKineticEnergyLoss());
  SetHistoryOutputValue("TotPressureLoss_Stage", TurboStagePerf->GetTotalPressureLoss());
}

void CFlowCompOutput::WriteTurboSpanwisePerformance(std::shared_ptr<CTurboOutput> TurboPerf, CGeometry *geometry, CConfig **config, unsigned short val_iZone) {

  string inMarker_Tag, outMarker_Tag, inMarkerTag_Mix;
  unsigned short nZone       = config[val_iZone]->GetnZone();

  unsigned short iDim, iSpan;

  unsigned long iExtIter = config[val_iZone]->GetOuterIter();
  const su2double* SpanWiseValuesIn, *SpanWiseValuesOut;
  ofstream file;
  string spanwise_performance_filename;

  auto BladePerformance = TurboPerf->GetBladesPerformances();

  /*--- Start of write file turboperformance spanwise ---*/
  SpanWiseValuesIn = geometry->GetSpanWiseValue(INFLOW);
  SpanWiseValuesOut = geometry->GetSpanWiseValue(OUTFLOW);

  /*--- Writing Span wise inflow thermodynamic quantities. ---*/
  spanwise_performance_filename = "TURBOMACHINERY/inflow_spanwise_thermodynamic_values.dat";
  if (nZone > 1) {
    spanwise_performance_filename.append("_" + std::to_string(val_iZone) + ".dat");
  } else {
    spanwise_performance_filename.append(".dat");
  }
  file.open (spanwise_performance_filename.data(), ios::out | ios::trunc);
  file.setf(ios::scientific);
  file.precision(12);

  file << "TITLE = \"Inflow Spanwise Thermodynamic Values. iOuterIter = " << iExtIter << " \"" << endl;
  file << "VARIABLES =" << endl;

  file.width(30); file << "\"SpanWise Value[m]\"";
  file.width(15); file << "\"iSpan\"";
  file.width(30); file << "\"Pressure[Pa]\"";
  file.width(30); file << "\"TotalPressure[Pa]\"";
  file.width(30); file << "\"Temperature[K]\"";
  file.width(30); file << "\"TotalTemperature[K]\"";
  file.width(30); file << "\"Enthalpy[J]\"";
  file.width(30); file << "\"TotalEnthalpy[J]\"";
  file.width(30); file << "\"Density[kg/m3]\"";
  file.width(30); file << "\"Entropy[J/K]\"";
  file << endl;

  for(iSpan = 0; iSpan < config[val_iZone]->GetnSpanWiseSections(); iSpan++){
    const auto& BladePerf = BladePerformance.at(val_iZone).at(iSpan);

    file.width(30); file << SpanWiseValuesIn[iSpan];
    file.width(15); file << iSpan;
    file.width(30); file << BladePerf->GetInletState().GetPressure()*config[ZONE_0]->GetPressure_Ref();
    file.width(30); file << BladePerf->GetInletState().GetTotalPressure()*config[ZONE_0]->GetPressure_Ref();
    file.width(30); file << BladePerf->GetInletState().GetTemperature()*config[ZONE_0]->GetTemperature_Ref();
    file.width(30); file << BladePerf->GetInletState().GetTotalTemperature()*config[ZONE_0]->GetTemperature_Ref();
    file.width(30); file << BladePerf->GetInletState().GetEnthalpy()*config[ZONE_0]->GetEnergy_Ref();
    file.width(30); file << BladePerf->GetInletState().GetTotalEnthalpy()*config[ZONE_0]->GetEnergy_Ref();
    file.width(30); file << BladePerf->GetInletState().GetDensity()*config[ZONE_0]->GetDensity_Ref();
    file.width(30); file << BladePerf->GetInletState().GetEntropy()*config[ZONE_0]->GetEnergy_Ref()/config[ZONE_0]->GetTemperature_Ref();
  }

  file.close();

  /*--- Writing Span wise outflow thermodynamic quantities. ---*/
  spanwise_performance_filename = "TURBOMACHINERY/outflow_spanwise_thermodynamic_values.dat";
  if (nZone > 1) {
    spanwise_performance_filename.append("_" + std::to_string(val_iZone) + ".dat");
  } else {
    spanwise_performance_filename.append(".dat");
  }
  file.open (spanwise_performance_filename.data(), ios::out | ios::trunc);
  file.setf(ios::scientific);
  file.precision(12);

  file << "TITLE = \"Outflow Span-wise Thermodynamic Values. iOuterIter = " << iExtIter << " \"" << endl;
  file << "VARIABLES =" << endl;

  file.width(30); file << "\"SpanWise Value[m]\"";
  file.width(15); file << "\"iSpan\"";
  file.width(30); file << "\"Pressure[Pa]\"";
  file.width(30); file << "\"TotalPressure[Pa]\"";
  file.width(30); file << "\"Temperature[K]\"";
  file.width(30); file << "\"TotalTemperature[K]\"";
  file.width(30); file << "\"Enthalpy[J]\"";
  file.width(30); file << "\"TotalEnthalpy[J]\"";
  file.width(30); file << "\"Density[kg/m3]\"";
  file.width(30); file << "\"Entropy[J/K]\"";
  file << endl;


  for(iSpan = 0; iSpan < config[val_iZone]->GetnSpanWiseSections(); iSpan++){
    const auto& BladePerf = BladePerformance.at(val_iZone).at(iSpan);

    file.width(30); file << SpanWiseValuesOut[iSpan];
    file.width(15); file << iSpan;
    file.width(30); file << BladePerf->GetOutletState().GetPressure()*config[ZONE_0]->GetPressure_Ref();
    file.width(30); file << BladePerf->GetOutletState().GetTotalPressure()*config[ZONE_0]->GetPressure_Ref();
    file.width(30); file << BladePerf->GetOutletState().GetTemperature()*config[ZONE_0]->GetTemperature_Ref();
    file.width(30); file << BladePerf->GetOutletState().GetTotalTemperature()*config[ZONE_0]->GetTemperature_Ref();
    file.width(30); file << BladePerf->GetOutletState().GetEnthalpy()*config[ZONE_0]->GetEnergy_Ref();
    file.width(30); file << BladePerf->GetOutletState().GetTotalEnthalpy()*config[ZONE_0]->GetEnergy_Ref();
    file.width(30); file << BladePerf->GetOutletState().GetDensity()*config[ZONE_0]->GetDensity_Ref();
    file.width(30); file << BladePerf->GetOutletState().GetEntropy()*config[ZONE_0]->GetEnergy_Ref()/config[ZONE_0]->GetTemperature_Ref();
  }

  file.close();

  /*--- Writing Span wise inflow kinematic quantities. ---*/
  spanwise_performance_filename = "TURBOMACHINERY/inflow_spanwise_kinematic_values.dat";
  if (nZone > 1) {
    spanwise_performance_filename.append("_" + std::to_string(val_iZone) + ".dat");
  } else {
    spanwise_performance_filename.append(".dat");
  }
  file.open (spanwise_performance_filename.data(), ios::out | ios::trunc);
  file.setf(ios::scientific);
  file.precision(12);

  file << "TITLE = \"Inflow Span-wise Kinematic Values. iOuterIter = " << iExtIter << " \"" << endl;
  file << "VARIABLES =" << endl;

  file.width(30); file << "\"SpanWise Value[m]\"";
  file.width(15); file << "\"iSpan\"";
  file.width(30); file << "\"Normal Mach[-]\"";
  file.width(30); file << "\"Tangential Mach[-]\"";
  if (geometry->GetnDim() == 3) {
    file.width(30); file << "\"3rd Component Mach[-]\"";
  };
  file.width(30); file << "\"Mach Module[-]\"";
  file.width(30); file << "\"Normal Velocity[m/s]\"";
  file.width(30); file << "\"Tangential Velocity[m/s]\"";
  if (geometry->GetnDim() == 3) {
    file.width(30); file << "\"3rd Component Velocity[m/s]\"";
  };
  file.width(30); file << "\"Velocity Module[m/s]\"";
  file.width(30); file << "\"Absolute Flow Angle[deg]\"";
  file.width(30); file << "\"Relative Flow Angle[deg]\"";
  file << endl;


  for(iSpan = 0; iSpan < config[val_iZone]->GetnSpanWiseSections(); iSpan++){
    const auto& BladePerf = BladePerformance.at(val_iZone).at(iSpan);

    file.width(30); file << SpanWiseValuesIn[iSpan];
    file.width(15); file << iSpan;
    for (iDim = 0; iDim < geometry->GetnDim(); iDim++){
      file.width(30); file << BladePerf->GetInletState().GetMach()[iDim];
    }
    file.width(30); file << BladePerf->GetInletState().GetMachValue();
    for (iDim = 0; iDim < geometry->GetnDim(); iDim++){
      file.width(30); file << BladePerf->GetInletState().GetVelocity()[iDim]*config[ZONE_0]->GetVelocity_Ref();
    }
    file.width(30); file << BladePerf->GetInletState().GetVelocityValue()*config[ZONE_0]->GetVelocity_Ref();
    // This captures NaNs
    if(isnan(BladePerf->GetInletState().GetAbsFlowAngle())){
      file.width(30); file << "0.0000";
    }
    else {
      file.width(30); file << BladePerf->GetInletState().GetAbsFlowAngle()*180.0/PI_NUMBER;
    }
    if(isnan(BladePerf->GetInletState().GetFlowAngle())){
      file.width(30); file << "0.0000";
    }
    else{
      file.width(30); file << BladePerf->GetInletState().GetFlowAngle()*180.0/PI_NUMBER;
    }
    file << endl;
  }

  file.close();

  /*--- Writing Span wise outflow thermodynamic quantities. ---*/
  spanwise_performance_filename = "TURBOMACHINERY/outflow_spanwise_kinematic_values";
  if (nZone > 1) {
    spanwise_performance_filename.append("_" + std::to_string(val_iZone) + ".dat");
  } else {
    spanwise_performance_filename.append(".dat");
  }
  file.open (spanwise_performance_filename.data(), ios::out | ios::trunc);
  file.setf(ios::scientific);
  file.precision(12);

  file << "TITLE = \"Outflow Span-wise Kinematic Values. iOuterIter = " << iExtIter << " \"" << endl;
  file << "VARIABLES =" << endl;

  file.width(30); file << "\"SpanWise Value[m]\"";
  file.width(15); file << "\"iSpan\"";
  file.width(30); file << "\"Normal Mach[-]\"";
  file.width(30); file << "\"Tangential Mach[-]\"";
  if (geometry->GetnDim() == 3) {
    file.width(30); file << "\"3rd Component Mach[-]\"";
  };
  file.width(30); file << "\"Mach Module[-]\"";
  file.width(30); file << "\"Normal Velocity[m/s]\"";
  file.width(30); file << "\"Tangential Velocity[m/s]\"";
  if (geometry->GetnDim() == 3) {
    file.width(30); file << "\"3rd Component Velocity[m/s]\"";
  };
  file.width(30); file << "\"Velocity Module[m/s]\"";
  file.width(30); file << "\"Absolute Flow Angle[deg]\"";
  file.width(30); file << "\"Relative Flow Angle[deg]\"";
  file << endl;


  for(iSpan = 0; iSpan < config[val_iZone]->GetnSpanWiseSections(); iSpan++){
    const auto& BladePerf = BladePerformance.at(val_iZone).at(iSpan);

    file.width(30); file << SpanWiseValuesOut[iSpan];
    file.width(15); file << iSpan;
    for (iDim = 0; iDim < geometry->GetnDim(); iDim++){
      file.width(30); file << BladePerf->GetOutletState().GetMach()[iDim];
    }
    file.width(30); file << BladePerf->GetInletState().GetMachValue();
    for (iDim = 0; iDim < geometry->GetnDim(); iDim++){
      file.width(30); file << BladePerf->GetOutletState().GetVelocity()[iDim]*config[ZONE_0]->GetVelocity_Ref();
    }
    file.width(30); file << BladePerf->GetInletState().GetVelocityValue()*config[ZONE_0]->GetVelocity_Ref();
    if(isnan(BladePerf->GetInletState().GetAbsFlowAngle())){
      file.width(30); file << "0.0000";
    }
    else {
      file.width(30); file << BladePerf->GetOutletState().GetAbsFlowAngle()*180.0/PI_NUMBER;
    }
    if(isnan(BladePerf->GetInletState().GetAbsFlowAngle())){
      file.width(30); file << "0.0000";
    }
    else{
      file.width(30); file << BladePerf->GetOutletState().GetFlowAngle()*180.0/PI_NUMBER;
    }
    file << endl;
  }

  file.close();
}
