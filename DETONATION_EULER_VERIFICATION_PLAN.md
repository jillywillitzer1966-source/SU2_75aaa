# DETONATION_EULER Verification Plan

## 1. Frozen Interface

This document freezes the current `DETONATION_EULER` Phase-1 MVP interface and
verification gates. It is intentionally verification-first:

- first verify chemistry energy closure in homogeneous 0D
- then verify homogeneous preservation
- then verify inert Euler transport
- only after these pass should reactive 1D / detonation claims be considered

Current limitations that must be stated in every report:

- `CDetonationNumerics` currently uses `DETONATION_LLF`, a local Lax-Friedrichs /
  Rusanov-like placeholder flux
- true AUSM+UP2 is **not** implemented yet
- no detonation validation exists yet

## 2. State Vector Ordering

The current conservative state vector is:

`U = [rho_0, rho_1, ..., rho_(Ns-1), rho*u, rho*v, rho*w, rhoE]`

For 2D, the stored order is:

`U = [rho_s..., rho*u, rho*v, rhoE]`

Definitions:

- `rho_s` : species partial density, `kg/m^3`
- `rho` : mixture density, `kg/m^3`, computed as `sum_s rho_s`
- `rho*u_i` : momentum density, `kg/(m^2 s)`
- `rhoE` : total energy density, `J/m^3`
- `E = e + 0.5 |u|^2`

No `rhoE_ve`, `Tve`, or other non-equilibrium modal energies are stored.

## 3. Primitive Variable Ordering

The current primitive/state-recovery vector is:

`V = [rho_0, rho_1, ..., rho_(Ns-1), T, u, v, w, p, rho, h, a]`

For 2D, the stored order is:

`V = [rho_s..., T, u, v, p, rho, h, a]`

Definitions:

- `T` : mixture temperature, `K`
- `u, v, w` : velocity components, `m/s`
- `p` : pressure, `Pa`
- `rho` : mixture density, `kg/m^3`
- `h` : mixture specific total enthalpy proxy stored from recovered conservative state, `J/kg`
- `a` : frozen-state sound speed from Cantera or fallback closure, `m/s`

## 4. Stored Quantity Meaning and Units

| quantity | meaning | unit |
| --- | --- | --- |
| `rho_s` | species partial density | `kg/m^3` |
| `rho` | mixture density, `sum rho_s` | `kg/m^3` |
| `u, v, w` | velocity components | `m/s` |
| `rhoE` | total energy density | `J/m^3` |
| `T` | mixture temperature | `K` |
| `p` | mixture pressure | `Pa` |
| `e` | mixture specific internal energy | `J/kg` |
| `h_s` | species specific enthalpy used by chemistry diagnostic | `J/kg` |
| `u_s` | species specific internal energy used by chemistry ODE | `J/kg` |
| `wdot_s` | species mass production rate after conversion from Cantera molar rates | `kg/(m^3 s)` |
| `qdot` | diagnostic volumetric heat-release proxy, `-sum(h_s * wdot_s)` | `W/m^3` |

## 5. Exact Energy Convention

- Cantera mixture internal energy is obtained from `ThermoPhase::intEnergy_mass()`, which is mass-based.
- Cantera species internal energies are obtained from `getPartialMolarIntEnergies()`, which returns molar values; the current solver divides by molecular weight to convert to `J/kg`.
- Cantera species enthalpies are obtained from `getPartialMolarEnthalpies()`, which returns molar values; the current solver divides by molecular weight to convert to `J/kg`.
- The Cantera thermodynamic state is tied to the mechanism reference-state thermochemistry. In the current implementation this means the reported internal energies / enthalpies are thermochemical quantities consistent with the mechanism thermo polynomials, including formation-energy offsets embedded in the reference-state model.
- Cantera net production rates are obtained from `Kinetics::getNetProductionRates()`, which returns molar production rates; the current solver multiplies by molecular weight to convert to `kg/(m^3 s)`.
- `qdot` is defined as `-sum(h_s * wdot_s)`.
- `qdot` is currently a diagnostic output only; the chemistry update is applied through the evolved `(T, Y_s)` state and then `rhoE = rho * (e(T,Y) + 0.5 |u|^2)` is reconstructed.

## 6. Cantera Unit Audit

| quantity | Cantera API source | raw Cantera unit | SU2 internal unit | conversion applied | verified yes/no |
| --- | --- | --- | --- | --- | --- |
| net production rates | `Kinetics::getNetProductionRates` | `kmol/(m^3 s)` | `kg/(m^3 s)` | multiply by molecular weight `kg/kmol` | yes |
| species enthalpies | `ThermoPhase::getPartialMolarEnthalpies` | `J/kmol` | `J/kg` | divide by molecular weight `kg/kmol` | yes |
| species internal energies | `ThermoPhase::getPartialMolarIntEnergies` | `J/kmol` | `J/kg` | divide by molecular weight `kg/kmol` | yes |
| mixture internal energy | `ThermoPhase::intEnergy_mass` | `J/kg` | `J/kg` | none | yes |
| mixture enthalpy | `ThermoPhase::enthalpy_mass` | `J/kg` | `J/kg` | none | yes |
| mixture `c_v` | `ThermoPhase::cv_mass` | `J/(kg K)` | `J/(kg K)` | none | yes |
| density | `ThermoPhase::density` | `kg/m^3` | `kg/m^3` | none | yes |
| pressure | `ThermoPhase::pressure` | `Pa` | `Pa` | none | yes |
| sound speed | `ThermoPhase::soundSpeed` | `m/s` | `m/s` | none | yes |
| `qdot` | computed in `CCanteraReactiveGas::ComputeChemistrySourceTerms` | `W/m^3` after conversion | `W/m^3` | `qdot = -sum(h_s * wdot_s)` | yes |

## 7. Current Verification Gates

### 7.1 0D homogeneous chemistry

Required checks:

- `CHEM_FAIL_CELLS = 0`
- `CHEM_E_REL_MAX < 1e-8` if achievable
- if `CHEM_E_REL_MAX > 1e-6`, mark the case failed
- all `Y_s >= 0` after clipping / renormalization
- `sum(Y_s) = 1` within `1e-12`
- compare cell-averaged SU2 outputs against an independent Cantera 0D reference

### 7.2 Homogeneous preservation

Required checks on the uniform 5x5 thermal bath:

- `max(T) - min(T)` near roundoff or explainable by MPI / CSV precision
- `max(p) - min(p)` near roundoff or explainable by MPI / CSV precision
- `max(rho) - min(rho)` near roundoff or explainable by MPI / CSV precision
- `max(Y_s) - min(Y_s)` near roundoff or explainable by MPI / CSV precision

### 7.3 Inert Euler verification

Required tests before reactive 1D claims:

- uniform free-stream preservation
- open supersonic inlet/outlet uniform-flow preservation
- 1D Sod-like inert Euler regression
- wall / symmetry reflection check

## 8. Current Verified Status

1. Build status
   Pass.
   `ninja -C /home/jmy/detonationFoam/SU2_deto/build install` completed successfully after the open-BC fix and verification-only left/right initialization addition.
2. Runtime smoke status
   Pass for the exercised `np=2` verification cases.
   Verified runs this round:
   - `detonation_homogeneous_drift_debug/chemistry_only.cfg`
   - `detonation_homogeneous_drift_debug/flow_only_zero_velocity.cfg`
   - `detonation_homogeneous_drift_debug/flow_plus_chemistry.cfg`
   - `detonation_open_bc_debug/supersonic_inlet_outlet_uniform_x.cfg`
   - `detonation_open_bc_debug/supersonic_outlet_only_uniform_x.cfg`
   - `detonation_open_bc_debug/closed_symmetry_uniform_flow.cfg`
   - `detonation_0d_h2o2n2/thermalbath_h2o2n2.cfg`
   - `detonation_inert_euler_tests/uniform_freestream.cfg`
   - `detonation_inert_euler_tests/symmetry_reflection.cfg`
   - `detonation_inert_euler_tests/sod_like_inert.cfg`
3. Chemistry backend status
   Pass at the unit-audit level.
   Cantera thermochemistry / kinetics backend remains wired into `DETONATION_EULER`, with verified unit conversions:
   - `wdot_s`: `kmol/(m^3 s)` -> `kg/(m^3 s)`
   - `h_s`, `u_s`: `J/kmol` -> `J/kg`
   - `c_v`: `J/(kg K)` unchanged
4. CVODE chemistry status
   Pass for the homogeneous verification window.
   - `CHEM_FAIL_CELLS = 0`
   - `CHEM_SUBSTEPS_MAX = 1`
5. Energy closure status
   Pass.
   Latest homogeneous H2/O2/N2 checks:
   - chemistry-only case: `CHEM_E_REL_MAX = 6.490752987e-11`, `CHEM_E_ABS_MAX = 1.210933551e-4 J/kg`
   - flow+chemistry case: `CHEM_E_REL_MAX = 6.490752987e-11`, `CHEM_E_ABS_MAX = 1.210933551e-4 J/kg`
   Interpretation:
   the chemistry substep energy closure remains trustworthy after the flow-path fixes.
6. Homogeneous 0D validation status
   Pass for the closed 5x5 thermal-bath target that motivated this debug round.
   Quantified SU2 vs independent Cantera 0D reference error:
   - `max |T_SU2 - T_ref| = 7.3132472051 K`
   - `max relative T error = 3.5656342540e-3`
   - `max relative pressure error = 3.4078197326e-3`
   - `max relative density error = 3.3084948197e-10`
   - `max relative internal-energy error = 4.4325498563e-11`
   - `max absolute species error = 1.2756860238e-2`
   Homogeneous-preservation before/after table:

   | case | metric | before fix | after fix |
   | --- | --- | --- | --- |
   | full 5x5 bath | `max(T)-min(T)` | `7.581387235e-1 K` | `4.320099833e-11 K` |
   | full 5x5 bath | `max(p)-min(p)` | `4.335437189e+1 Pa` | `2.124579623e-9 Pa` |
   | full 5x5 bath | `max(rho)-min(rho)` | `9.750648171e-5 kg/m^3` | `1.387778781e-16 kg/m^3` |

   Closed-domain verification matrix after the fix:
   - `CHEMISTRY_ONLY`:
     - `DetSpanTemperature_max_over_time = 3.251443559e-11 K`
     - `DetSpanPressure_max_over_time = 1.571606845e-9 Pa`
     - `DetSpanDensity_max_over_time = 1.387778781e-16 kg/m^3`
   - `FLOW_ONLY_ZERO_VELOCITY`:
     - `DetSpanTemperature_max_over_time = 0`
     - `DetSpanPressure_max_over_time = 0`
     - `DetSpanDensity_max_over_time = 0`
   - `FLOW_PLUS_CHEMISTRY`:
     - `DetSpanTemperature_max_over_time = 4.320099833e-11 K`
     - `DetSpanPressure_max_over_time = 2.124579623e-9 Pa`
     - `DetSpanDensity_max_over_time = 1.387778781e-16 kg/m^3`
7. Inert Euler validation status
   Partial, but materially improved.
   Verified pass cases:
   - inert uniform free-stream:
     - `DetSpanTemperature_max_over_time = 0`
     - `DetSpanPressure_max_over_time = 0`
     - `DetSpanDensity_max_over_time = 0`
   - open supersonic inlet/outlet uniform flow:
     - before fix:
       - `DetSpanTemperature_max_over_time = 2.723629731 K`
       - `DetSpanPressure_max_over_time = 9.199051683e+2 Pa`
       - `DetSpanDensity_max_over_time = 9.824628222e-9 kg/m^3`
     - after fix:
       - `DetSpanTemperature_max_over_time = 0`
       - `DetSpanPressure_max_over_time = 0`
       - `DetSpanDensity_max_over_time = 0`
       - total mass drift = `0`
       - total x-momentum drift = `0`
       - total energy drift = `0`
   - outlet-only uniform flow:
     - all scalar spans = `0`
     - total mass / x-momentum / energy drifts = `0`
   - closed symmetry uniform flow:
     - all scalar spans = `0`
     - `max(rho*u)-min(rho*u) = 5.456968211e-16`
   - symmetry reflection runtime path remains stable

   `farfield_uniform_flow` status:
   - not a drift failure
   - currently blocked by the SU2 framework/config restriction:
     `Riemann Boundary conditions or Giles must be used outlet with Not Ideal Compressible Fluids`

   1D inert Sod-like status:
   - no longer blocked
   - verification-only left/right initialization path added through `DETONATION_LEFT_RIGHT_INIT = YES`
   - current results:
     - `rho_l1_error = 1.7473751996717118e-02`
     - `u_l1_error = 2.0434989048805368`
     - `p_l1_error = 1.220781857937722e+03 Pa`
     - observed left-going disturbance = `true`
     - observed right-going disturbance = `true`
   Interpretation:
   - the current inert Euler path is good enough for indexing and wave-direction sanity checks
   - it is not yet a final closed-domain conservation benchmark because the current Sod-like setup uses open `x` boundaries
8. Known failures
   - The current convective path is `DETONATION_LLF`, not AUSM+UP2.
   - `DETONATION_EULER` still cannot claim 1D reactive shock / detonation validation.
   - `farfield_uniform_flow` remains framework-blocked, so the FAR_FIELD path is not yet exercised as a pass case.
   - the current Sod-like inert case uses a verification-only initialization path and open `x` boundaries
   - restart-based conservation summaries should still be interpreted carefully for `DETONATION_EULER`
9. Next recommended task
   Do not proceed to reactive 1D yet.
   If verification continues, the next worthwhile step is to strengthen the inert Sod-like benchmark into a more conservative and less boundary-sensitive reference case, then revisit higher-accuracy inert Euler validation before any reactive 1D work.

## 9. 0D homogeneous chemistry verification milestone

Accepted milestone values:

- `CHEM_FAIL_CELLS = 0`
- `CHEM_SUBSTEPS_MAX = 1`
- `CHEM_E_REL_MAX = 6.490752987e-11`
- closed homogeneous bath spans:
  - `T span = 4.320099833e-11 K`
  - `p span = 2.124579623e-9 Pa`
  - `rho span = 1.387778781e-16 kg/m^3`
- open uniform-flow preservation: PASS

Interpretation:

- this milestone proves 0D chemistry closure and homogeneous-preservation
  integrity
- it does **not** prove detonation validation

## 10. Inert Euler benchmark milestone (2026-05-11)

### 10.1 Benchmark directory

- `/home/jmy/detonationFoam/SU2_test_case/detonation_inert_euler_benchmark`

### 10.2 Native output-chain status

Tested formats:

- `RESTART_ASCII`: success
- `PARAVIEW_ASCII`: success
- `PARAVIEW` XML: success

Audited fields:

- `Density`
- `Momentum`
- `Energy`
- `RhoS_i`
- `MassFrac_i`
- `Pressure`
- `Temperature`
- `Mach`
- `Sound_Speed`

Medium-grid audit highlights:

- max `|Density_vtk - Density_restart| = 4.855708919393464e-07`
- max `|Momentum_x_vtk - Momentum_x_restart| = 4.973616628944910e-05`
- max `|Energy_vtk - Energy_restart| = 4.999768556444906e-03`
- no missing required scalar/vector fields

Conclusion:

- native ParaView-readable volume output is working for `DETONATION_EULER`
- output indexing is consistent with the detonation state layout

### 10.3 Closed Sod-like inert benchmark

Grids:

- `100`
- `200`
- `400`

Observed wave directions:

- shock: right
- contact: right
- rarefaction: left

Approximate positions on the fine grid:

- observed shock position: `0.9425`
- exact shock position: `0.9421560736970935`
- observed contact position: `0.735`
- exact contact position: `0.734042443559385`

L1 error trend:

- `rho`: `2.896291609688977e-02 -> 1.897256640696782e-02 -> 1.186461481826517e-02`
- `u`: `1.554952596849264e+01 -> 8.785950451027794e+00 -> 4.769993969811145e+00`
- `p`: `2.727353012866930e+03 -> 1.527654416907767e+03 -> 8.670998496924260e+02`

Grid-refinement conclusion:

- `rho_error_decreases = true`
- `u_error_decreases = true`
- `p_error_decreases = true`

### 10.4 Closed-domain conservation status

Restart-based integrated drifts:

- `sod_closed_100`
  - `max_rel_mass_drift = 6.029637562376491e-14`
  - `max_rel_energy_drift = 6.992447226406705e-09`
- `sod_closed_200`
  - `max_rel_mass_drift = 1.130013175676776e-14`
  - `max_rel_energy_drift = 6.348522235641731e-09`
- `sod_closed_400`
  - `max_rel_mass_drift = 1.359553766349214e-14`
  - `max_rel_energy_drift = 6.074233513421270e-09`

Interpretation:

- total mass conservation is very strong
- total energy drift is small and consistent with output-derived integration
  from nodal restart fields, not with a physics inconsistency
- total x-momentum is not used as a strict conservation metric because the
  closed Euler-wall tube experiences wall pressure forces

### 10.5 Verified-case package

Created and rerun successfully:

- `/home/jmy/detonationFoam/SU2_test_case/verified_detonation_euler_cases`

Packaged cases:

- `00_chemistry_only_0d`
- `01_closed_flow_plus_chemistry_0d`
- `02_open_uniform_supersonic_flow`
- `03_inert_uniform_freestream`
- `04_inert_symmetry_reflection`
- `05_inert_sod_like_benchmark`

Package root runner:

```bash
cd /home/jmy/detonationFoam/SU2_test_case/verified_detonation_euler_cases
./run_all_verified.sh
```

### 10.6 Updated status

Current conclusion for the inert Euler stage:

- **A. PASS**
  inert Euler path and output chain are sufficiently verified; the verified
  case package is created; the code may proceed to reactive 1D validation in a
  later round.

This does **not** mean:

- reactive 1D has already been validated
- detonation has been validated

## 11. Inert Euler and output-chain verification milestone

Accepted milestone values:

- native output-chain audit:
  - `RESTART_ASCII`: success
  - `PARAVIEW_ASCII`: success
  - `PARAVIEW XML (.vtu)`: success
- closed Sod conservation:
  - max relative mass drift on `100/200/400` grids:
    `6.03e-14 / 1.13e-14 / 1.36e-14`
  - max relative energy drift on `100/200/400` grids:
    `6.99e-09 / 6.35e-09 / 6.07e-09`
- Sod exact-comparison trend:
  - `rho_L1: 2.896e-02 -> 1.897e-02 -> 1.186e-02`
  - `u_L1: 1.555e+01 -> 8.786e+00 -> 4.770e+00`
  - `p_L1: 2.727e+03 -> 1.528e+03 -> 8.671e+02`
- `verified_detonation_euler_cases`:
  - all `6` cases `PASS`

Interpretation:

- this milestone proves the inert Euler path and native output chain are
  trustworthy enough to start reactive 1D verification
- this is **not** reactive 1D validation
- this is **not** detonation validation
- the current flux is still `DETONATION_LLF`

## 12. Reactive 1D validation milestone (2026-05-11)

Preflight status:

- `/home/jmy/detonationFoam/SU2_test_case/detonation_reactive_1d_validation/run_preflight_verified.sh`
- verified suite result: `PASS`

Reactive uniform 1D regression:

- `CHEM_FAIL_CELLS = 0`
- `CHEM_SUBSTEPS_MAX = 1`
- `CHEM_E_REL_MAX = 6.490865307e-11`
- spatial spans remained near roundoff
- SU2 vs independent Cantera 0D reference stayed consistent with the accepted
  homogeneous 0D milestone

Reactive shock-tube baseline:

- final adopted setup:
  - same `H2/O2/N2` mixture as the verified homogeneous chemistry cases
  - driven gas: `P = 101325 Pa`, `T = 900 K`
  - driver gas: `P = 5.06625e7 Pa`, `T = 900 K`
- frozen comparison shows clear shock-triggered chemistry
  - `peak_temperature_reactive = 3566.161 K`
  - `peak_temperature_frozen = 1923.093 K`
  - `h2_consumption_max = 2.75474518e-2`
  - `h2o_production_max = 2.432286e-1`
- native reactive output-chain audit:
  - `RESTART_ASCII`: success
  - `PARAVIEW_ASCII`: success
  - `PARAVIEW XML (.vtu)`: success

Current blockers that prevent a full reactive-1D pass:

- strong reactive cells reach only
  `CHEM_E_REL_MAX = 3.633713259e-07 -> 5.170493412e-07 -> 5.999206678e-07`
  across `dt_base -> dt_half -> dt_quarter`
- Lie-splitting time-step sensitivity is significant
  - `significant = true`
  - current recommendation:
    `Strang splitting should be considered before HLLC-P.`

Current conclusion:

- **B. PARTIAL PASS**
  reactive 1D coupling and reactive output chain are working and physically
  informative, but the strong-shock chemistry-energy diagnostic and time-step
  sensitivity are not yet clean enough for detonation validation

## 13. Reactive 1D cleanup milestone (2026-05-11)

Starting point accepted for this round:

- reactive uniform 1D: pass
- reactive shock tube: usable baseline
- strong reactive cells:
  - `dt_base CHEM_E_REL_MAX = 3.633713259e-07`
  - `dt_half CHEM_E_REL_MAX = 5.170493412e-07`
  - `dt_quarter CHEM_E_REL_MAX = 5.999206678e-07`
- Lie time-step sensitivity:
  - `base / half / quarter wave speed = 852.816 / 863.903 / 863.770 m/s`
  - `significant = true`

Cleanup-round diagnostic findings:

- worst-cell chemistry-energy audit:
  - `CHEM_E_REL_MAX = 3.633713259207901e-07`
  - `CHEM_E_REL_P99 = 3.1235375665248864e-12`
  - `CHEM_E_REL_P999 = 2.1707217956715169e-10`
  - `num_cells_gt_1e-7 = 3`
  - `num_cells_gt_1e-6 = 0`
  - no meaningful clipping or renormalization correction
- CVODE tolerance sweep:
  - tighter tolerances improve `P99/P999`
  - tighter tolerances do **not** materially reduce the worst-cell
    `CHEM_E_REL_MAX`
  - larger `CHEMISTRY_MAX_SUBSTEPS` does not help because
    `CHEM_SUBSTEPS_MAX` stays `1`

Optional splitting work completed:

- new config option:
  - `CHEMISTRY_SPLITTING = LIE / STRANG`
  - default remains `LIE`
- an initial Strang regression was traced to `Solution_Old` not being
  synchronized before the second flow half-step
- fix applied:
  - synchronize `Solution_Old <- Solution` before each Strang flow half-step
- after the fix:
  - verified suite re-run: `PASS`
  - uniform reactive 1D with `STRANG`: matches `LIE`
    - `CHEM_E_REL_MAX = 6.4908278669804706e-11`
    - `T span = 3.956301953e-11 K`

Current Strang shock-tube status:

- `strang_base`:
  - physically reasonable
  - `reaction_triggered_behind_shock = true`
  - `wave_speed = 864.1691916873255 m/s`
- `strang_half`:
  - `reaction_triggered_behind_shock = true`
  - `wave_speed = 863.9028351930965 m/s`
  - `h2_consumption_max = 0.0275375764`
- `strang_quarter`:
  - `reaction_triggered_behind_shock = true`
  - `wave_speed = 866.043374215041 m/s`
  - `h2_consumption_max = 0.0275346155`

Current conclusion:

- **B. PARTIAL PASS**
  the strong-cell chemistry-energy diagnostic is now understood well enough to
  rule out a widespread energy-closure bug, and the homogeneous `STRANG` path
  is repaired. `STRANG` also preserves a physically meaningful reactive
  shock-tube path across `dt_base / dt_half / dt_quarter`, with a reduced
  wave-speed spread relative to `LIE`, but its worst-cell chemistry
  diagnostics remain comparable to `LIE`.

Current recommendation:

- continue reactive 1D cleanup
- do **not** proceed to planar detonation, CJ/ZND validation, `HLLC-P`, or
  `AUSM+UP2`

## Reactive 1D cleanup milestone: B. PARTIAL PASS

- STRANG implemented as optional experimental splitting.
- `Solution_Old` synchronization before the second flow half-step was fixed.
- Verified suite PASS reconfirmed.
- `uniform_lie` / `uniform_strang` PASS reconfirmed.
- Reactive shock-tube STRANG cases completed successfully.
- STRANG wave-speed consistency better than LIE:
  - LIE spread = `11.086930 m/s`
  - STRANG spread = `2.140539 m/s`
- Strong-cell chemistry diagnostics remain sparse but non-zero:
  - `CHEM_E_REL_MAX = O(1e-7)`
  - `CHEM_E_REL_P99 = 3.1235375665e-12`
  - `CHEM_E_REL_P999 = 2.1707217957e-10`
  - `num_cells_gt_1e-7 = 3`
  - `num_cells_gt_1e-6 = 0`
- Accepted handoff:
  - proceed to a gated 1D planar detonation candidate with CJ reference
  - do not proceed to 2D / HLLC-P / AUSM+UP2 / NS / AMR / DLB

## 1D planar detonation candidate status (current round)

- New workspace:
  `/home/jmy/detonationFoam/SU2_test_case/detonation_1d_planar_candidate`
- CJ reference tool:
  - local C++ Cantera implementation
  - method: equilibrium Hugoniot + Rayleigh minimum-speed search
- Stoichiometric `H2/O2/N2`, `900 K / 101325 Pa` CJ reference:
  - `D_CJ = 1904.729230971 m/s`
  - `p_CJ = 5.253040514e+05 Pa`
  - `T_CJ = 3056.545645237 K`
  - `rho_CJ = 4.821462662e-01 kg/m^3`
- Quantified coarse STRANG pilot:
  - quasi-1D mesh `(701, 3, 1)`
  - `dt = 2.0e-8 s`
  - chemistry on
  - `DETONATION_LLF`
  - `CHEMISTRY_SPLITTING = STRANG`
- Front-speed measurement over the selected steady window:
  - pressure-gradient = `1730.0 m/s`
  - qdot-peak = `1730.0 m/s`
  - temperature-threshold = `1700.0 m/s`
  - detector spread = `30.0 m/s`
- CJ comparison on the coarse STRANG pilot:
  - relative speed error = `9.173441985%`
  - post-front pressure error = `13.741098780%`
  - post-front temperature error = `14.647209537%`
  - post-front density error = `28.207837186%`
- Chemistry / conservation diagnostics on the coarse pilot:
  - `CHEM_FAIL_CELLS = 0`
  - `CHEM_SUBSTEPS_MAX = 1`
  - `CHEM_E_REL_MAX = 5.358232982e-10`
  - `max_rel_mass_drift = 1.0465750198e-08`
  - `max_rel_energy_drift = 7.9838154518e-08`
- Native output-chain audit on the planar candidate coarse pilot:
  - `RESTART_ASCII`: success
  - `PARAVIEW_ASCII`: success
  - `PARAVIEW XML (.vtu)`: success
- Status:
  - candidate propagation is promising and now quantified against a CJ reference
  - full LIE vs STRANG baseline, dt sweep, and grid sweep are still incomplete
  - no detonation-validation claim may be made yet
  - current conclusion remains **B. PARTIAL PASS**

## 1D planar detonation candidate progress update

- New reproducible sweep entry points have been added:
  - `cases/03_planar_strang_dt_sweep/run.sh`
  - `cases/03_planar_strang_dt_sweep/postprocess.sh`
  - `cases/04_planar_strang_grid_sweep/run.sh`
  - `cases/04_planar_strang_grid_sweep/postprocess.sh`
- Coarse planar LIE vs STRANG baseline is now quantified:
  - LIE:
    - mean speed = `1723.333333333 m/s`
    - detector spread = `20.0 m/s`
    - CJ speed error = `9.173441985%`
    - `CHEM_E_REL_MAX = 6.015580761e-10`
  - STRANG:
    - mean speed = `1720.0 m/s`
    - detector spread = `30.0 m/s`
    - CJ speed error = `9.173441985%`
    - `CHEM_E_REL_MAX = 5.358232982e-10`
  - interpretation:
    the current LLF coarse-planar LIE and STRANG pilots are nearly identical in
    mean speed and CJ error; this comparison does not justify changing the
    default splitting mode.
- STRANG dt sweep progress:
  - `dt_base` complete:
    - mean speed = `1719.9999999999993 m/s`
    - pressure-gradient speed = `1729.9999999999993 m/s`
    - relative CJ speed error = `9.173441984803962%`
    - `CHEM_E_REL_MAX = 5.358232982e-10`
    - `CHEM_FAIL_CELLS = 0`
    - `max_rel_mass_drift = 1.2854671287523957e-08`
    - `max_rel_energy_drift = 7.983815451827048e-08`
  - `dt_half` complete:
    - mean speed = `1713.333333333333 m/s`
    - pressure-gradient speed = `1710.0000000000002 m/s`
    - relative CJ speed error = `10.223459996540285%`
    - `CHEM_E_REL_MAX = 2.281611268e-07`
    - `CHEM_FAIL_CELLS = 0`
    - `max_rel_mass_drift = 1.0693266435527115e-08`
    - `max_rel_energy_drift = 9.200600155264638e-08`
- Current interpretation:
  - the planar candidate continues to produce a right-running detonation-like
    wave with quantified CJ error in the `9-10%` band for the current LLF
    baseline
  - dt sensitivity is still active and must not be dismissed
  - no chemistry failures have appeared, but `dt_half` shows that planar
    strong-cell diagnostics can be notably harsher than the original coarse
    pilot
- Status remains **B. PARTIAL PASS**:
  - coarse-planar baseline is now better organized and reproducible
  - `dt_base` and `dt_half` are quantified
  - `dt_quarter` and the dedicated planar grid sweep remain pending

## Planar detonation candidate sweep milestone: current partial status

Execution gate:

- `verified_detonation_euler_cases/run_all_verified.sh`: `PASS`
- because local CPU resources are limited, main planar-candidate cases should
  be run **one at a time** and monitored to completion

Completed STRANG dt sweep:

- `dt_base`
  - mean speed = `1719.9999999999993 m/s`
  - CJ speed error = `9.173441984803962%`
  - `CHEM_E_REL_MAX = 5.358232982e-10`
  - `CHEM_FAIL_CELLS = 0`
- `dt_half`
  - mean speed = `1713.333333333333 m/s`
  - CJ speed error = `10.223459996540285%`
  - `CHEM_E_REL_MAX = 2.281611267539021e-07`
  - `CHEM_E_REL_P99 = 3.5949496150092683e-13`
  - `CHEM_E_REL_P999 = 5.2144228045831814e-11`
  - `num_cells_gt_1e-7 = 3`
  - `num_cells_gt_1e-6 = 0`
  - `CHEM_FAIL_CELLS = 0`
- `dt_quarter`
  - mean speed = `1716.666666666666 m/s`
  - CJ speed error = `10.223459996540321%`
  - `CHEM_E_REL_MAX = 2.0226547600352885e-10`
  - `CHEM_E_REL_P99 = 4.3256043487773676e-13`
  - `CHEM_E_REL_P999 = 1.8792253500890592e-11`
  - `num_cells_gt_1e-7 = 0`
  - `num_cells_gt_1e-6 = 0`
  - `CHEM_FAIL_CELLS = 0`

dt-sweep summary:

- mean-speed spread = `6.666666666666288 m/s`
- relative mean-speed spread = `0.3883495145630849%`
- conclusion:
  `dt sensitivity is reasonably controlled for the current LLF planar candidate baseline.`

Completed STRANG grid sweep:

- `grid_401`
  - mean speed = `1715.0 m/s`
  - CJ speed error = `9.960955493606206%`
  - `CHEM_E_REL_MAX = 2.3180427185757515e-07`
  - `CHEM_E_REL_P99 = 8.093149206873554e-13`
  - `CHEM_E_REL_P999 = 4.6688586431918606e-11`
  - `num_cells_gt_1e-7 = 6`
  - `num_cells_gt_1e-6 = 0`
  - `CHEM_FAIL_CELLS = 0`
- `grid_701`
  - mean speed = `1713.333333333333 m/s`
  - CJ speed error = `10.223459996540285%`
  - `CHEM_E_REL_MAX = 2.281611267539021e-07`
  - `CHEM_E_REL_P99 = 3.5949496150092683e-13`
  - `CHEM_E_REL_P999 = 5.2144228045831814e-11`
  - `num_cells_gt_1e-7 = 3`
  - `num_cells_gt_1e-6 = 0`
  - `CHEM_FAIL_CELLS = 0`
- `grid_1401`
  - mean speed = `1736.6666666666667 m/s`
  - CJ speed error = `8.648432978935718%`
  - `CHEM_E_REL_MAX = 2.2650718278089985e-07`
  - `CHEM_E_REL_P99 = 2.925314665008863e-13`
  - `CHEM_E_REL_P999 = 5.59817549578758e-11`
  - `num_cells_gt_1e-7 = 12`
  - `num_cells_gt_1e-6 = 0`
  - `CHEM_FAIL_CELLS = 0`

grid-sweep summary:

- `cj_speed_error_monotone_improvement = false`
- conclusion:
  `CJ speed error does not improve monotonically with the current spatial refinement sweep; continue planar cleanup.`

Accepted conclusion after the completed local sweep:

- status remains **B. PARTIAL PASS**
- current LLF planar candidate shows a stable right-running wave with
  quantifiable CJ comparison
- time-step sensitivity is now reasonably controlled for this baseline
- the remaining blocker is the non-monotone grid/CJ trend, not chemistry
  failures
- do not claim detonation validation
- do not proceed to 2D / HLLC-P / AUSM+UP2 / NS / diffusion / AMR / DLB from
  this milestone alone

## LLF planar detonation candidate milestone: B. PARTIAL PASS

Accepted LLF baseline facts:

- verified suite: `PASS`
- the 1D planar candidate propagates stably with `DETONATION_LLF`
- `STRANG dt_base / dt_half / dt_quarter` are complete
- `401 / 701 / 1401` grid sweep is complete
- `CHEM_FAIL_CELLS = 0`
- `CHEM_E_REL_P99 / P999` are small
- strong-cell `CHEM_E_REL_MAX` spikes are sparse, not widespread

Key metrics:

- `dt_quarter`
  - speed mean = `1716.666666666666 m/s`
  - CJ speed error = `10.2234599965%`
  - `CHEM_E_REL_MAX = 2.0226547600e-10`
  - `CHEM_FAIL_CELLS = 0`
- `grid_1401`
  - speed mean = `1736.6666666666667 m/s`
  - CJ speed error = `8.6484329789%`
  - `CHEM_E_REL_MAX = 2.2650718278e-07`
  - `CHEM_FAIL_CELLS = 0`
- dt sweep:
  - mean-speed spread = `6.666666666666288 m/s`
  - relative spread = `0.3883495145630849%`
  - time-step sensitivity is basically controlled
- grid sweep:
  - `401 / 701 / 1401` complete
  - CJ error trend is not monotone

Conclusion:

- `DETONATION_LLF` remains the verified baseline/debug flux
- LLF is usable but too dissipative / not clean enough for `A. PASS`
- next verification stage may study 1D flux improvement with
  `DETONATION_HLLC`
- this is still not final detonation validation

## DETONATION_HLLC flux-improvement gate: current status

Purpose:

- study whether a less dissipative 1D Euler flux improves inert shock/contact
  resolution and planar detonation-candidate CJ-speed comparison
- retain `DETONATION_LLF` as baseline/debug flux
- do not change chemistry, splitting, energy convention, state-vector layout,
  Cantera/CVODE backend, or boundary physics

Implementation:

- `CONV_NUM_METHOD_FLOW= DETONATION_HLLC` is config-visible
- HLLC uses the current detonation conservative layout:
  `U = [rho_s..., rho*u_i, rhoE]`
- species are transported through the HLLC contact wave with side-specific mass
  fractions
- HLLC uses pressure/density/composition positivity fallback to LLF
- history diagnostics added:
  - `HLLC_FALLBACK_COUNT`
  - `HLLC_NEGATIVE_STATE_COUNT`
  - `HLLC_MIN_PRESSURE`
  - `HLLC_MIN_DENSITY`

Formula/unit check:

- report:
  `/home/jmy/detonationFoam/SU2_test_case/detonation_flux_improvement_1d/results/HLLC_FLUX_UNIT_CHECK.md`
- result: `PASS`
- uniform left/right states recover analytic physical flux for `+x`, `-x`,
  `+y`, and `-y` normals
- moving/stationary contact checks pass
- Sod left/right finite-flux check passes

Runtime finding:

- first HLLC open-uniform pilot preserved the uniform field but reported:
  - `HLLC_FALLBACK_COUNT = 52`
  - `HLLC_NEGATIVE_STATE_COUNT = 52`
- root cause:
  a positivity check assumed `rhoE > kinetic energy`
- this assumption is invalid for Cantera formation-energy thermochemistry,
  where physical mixture internal energy can be negative
- fix:
  remove the invalid `rhoE - kinetic` rejection and keep pressure, density,
  and composition positivity checks

Current required gate before further HLLC validation:

- rerun corrected HLLC open-uniform case with `np=2`
- rerun HLLC closed homogeneous flow+chemistry with `np=2`
- rerun HLLC inert uniform freestream with `np=2`
- rerun full verified suite after the positivity correction

Current conclusion:

- **B. PARTIAL PASS**
- corrected source is built and installed
- corrected MPI runtime verification is pending because the current tool
  session was blocked from launching additional `mpirun -np 2` jobs
- do not proceed to HLLC inert Sod, reactive shock tube, planar candidate,
  HLLC-P, or 2D until these corrected gates pass

## HLLC uniform gates: A. PASS

The corrected HLLC uniform gates were rerun with `np=2`, one main case at a
time.

### open_uniform_hllc

- max `T/p/rho` span = `0 / 0 / 0`
- max `rho*u/rho*v/rhoE` span = `0 / 7.275957614e-16 / 0`
- total mass drift = `0`
- total x-momentum drift = `0`
- total energy drift = `0`
- `HLLC_FALLBACK_COUNT = 0`
- `HLLC_NEGATIVE_STATE_COUNT = 0`
- `HLLC_MIN_PRESSURE = 101325 Pa`
- `HLLC_MIN_DENSITY = 0.8494838549 kg/m^3`

Result: `PASS`

### closed_flow_plus_chemistry_hllc

- max `T` span = `3.365130397e-11 K`
- max `p` span = `1.702574082e-09 Pa`
- max `rho` span = `2.220446049e-16 kg/m^3`
- `CHEM_FAIL_CELLS = 0`
- `CHEM_SUBSTEPS_MAX = 1`
- `CHEM_E_REL_MAX = 6.490740507e-11`
- `HLLC_FALLBACK_COUNT = 0`
- `HLLC_NEGATIVE_STATE_COUNT = 0`
- `HLLC_MIN_PRESSURE = 101315.5307 Pa`
- `HLLC_MIN_DENSITY = 0.1274225782 kg/m^3`

Result: `PASS`

### inert_uniform_hllc

- max `T/p/rho` span = `0 / 0 / 0`
- max `rho*u/rho*v/rhoE` span = `0 / 7.275957614e-16 / 0`
- total mass drift = `0`
- total x-momentum drift = `0`
- total energy drift = `0`
- `HLLC_FALLBACK_COUNT = 0`
- `HLLC_NEGATIVE_STATE_COUNT = 0`
- `HLLC_MIN_PRESSURE = 101325 Pa`
- `HLLC_MIN_DENSITY = 0.8494838549 kg/m^3`

Result: `PASS`

Uniform-gate conclusion:

- **A. PASS**
- the HLLC formation-energy positivity fix is verified at runtime
- HLLC uniform preservation is clean
- HLLC can proceed to smoke/manual inert Sod, reactive shock, and planar
  candidate preparation
- this does not validate HLLC inert Sod, reactive shock, planar detonation,
  HLLC-P, or 2D readiness

## HLLC smoke/manual package

Created:

`/home/jmy/detonationFoam/SU2_test_case/detonation_flux_improvement_1d/hllc_smoke_and_manual_package`

Smoke cases were run one at a time with `np=2`:

- `00_hllc_inert_sod_smoke`: `PASS`
- `01_hllc_reactive_shock_smoke`: `PASS`
- `02_hllc_planar_smoke`: `PASS`

Smoke diagnostics:

- all smoke cases:
  - `HLLC_FALLBACK_COUNT_max = 0`
  - `HLLC_NEGATIVE_STATE_COUNT_max = 0`
  - restart/VTK/VTU outputs generated
- reactive smoke cases:
  - `CHEM_FAIL_CELLS_max = 0`
  - reactive shock smoke `CHEM_E_REL_MAX = 1.871802097e-14`
  - planar smoke `CHEM_E_REL_MAX = 5.472412814e-12`

Smoke result:

- `PASS` for startup/diagnostics/output plumbing
- not a validation of full HLLC inert/reactive/planar behavior

## HLLC manual full-case verification status

The user manually completed the HLLC full cases in:

`/home/jmy/detonationFoam/SU2_test_case/detonation_flux_improvement_1d/hllc_smoke_and_manual_package/manual_cases`

The completed results were inspected offline. No new SU2 main runs were
launched for this closeout.

### Runtime and output chain

- `00_hllc_inert_sod_full`: solver exit success, restart/VTK/VTU present
- `01_hllc_reactive_shock_full`: solver exit success, restart/VTK/VTU present
- `02_hllc_planar_full`: solver exit success, restart/VTK/VTU present
- all cases: `HLLC_FALLBACK_COUNT_max = 0`
- all cases: `HLLC_NEGATIVE_STATE_COUNT_max = 0`
- reactive cases: `CHEM_FAIL_CELLS_max = 0`

Runtime/output conclusion: `PASS`.

### Inert Sod HLLC full gate

- max relative mixture-mass drift = `1.6954470568438713e-13`
- max relative total-energy drift = `7.65986955033685e-13`
- `rho_L1 = 2.064834693351159e-02`
- `u_L1 = 1.1271583334369236e+01`
- `p_L1 = 1.9815549383402727e+03`
- previous LLF 100-cell benchmark:
  `rho_L1 / u_L1 / p_L1 =
  2.896291609688977e-02 / 1.5549525968492645e+01 / 2.72735301286693e+03`
- HLLC preserves positivity and improves all three L1 errors relative to LLF

Gate result: `PASS`.

### Reactive shock HLLC full gate

- max relative mixture-mass drift = `3.1222077770829243e-07`
- max relative total-energy drift = `5.1282958947458685e-08`
- max `|sum(Y)-1| = 7.134900004146516e-08`
- `CHEM_E_REL_MAX_max = 1.372358265e-08`
- reactive peak pressure / temperature =
  `50220700 Pa / 3650.829 K`
- frozen peak pressure / temperature =
  `50076100 Pa / 1923.093 K`
- max H2 consumption = `2.79889259e-02`
- max O2 consumption = `2.23109672e-01`
- max H2O production = `2.488593e-01`
- peak `QDOT = 7.76327925e+13`
- the current conservative detector reports
  `reaction_triggered_behind_shock = false`
- the pressure-gradient shock detector is near `x = 0.0305`, while the main
  QDOT/reaction peak is near `x = 0.0315`

Gate result: `PARTIAL PASS`. The full run, output, chemistry, and HLLC
diagnostics are clean, but the shock/reaction ordering should be rechecked with
a refined detector or controlled higher-resolution/denser-output run.

### Planar HLLC full candidate

- max relative mixture-mass drift = `1.4674801884028627e-08`
- max relative total-energy drift = `1.6864209138150115e-07`
- max `|sum(Y)-1| = 6.828999987718021e-08`
- `CHEM_E_REL_MAX_max = 1.846877537e-10`
- pressure-gradient speed = `1750.869994559564 m/s`
- QDOT-peak speed = `1770.8799945545575 m/s`
- temperature-threshold speed = `1730.8549895620672 m/s`
- mean speed = `1750.8683262253962 m/s`
- detector spread = `40.02500499249027 m/s`
- CJ speed = `1904.729230970699 m/s`
- pressure-gradient CJ speed error = `8.077748475184812%`
- previous LLF medium-grid CJ speed error = `10.2234599965%`

Candidate result: `PARTIAL PASS`. HLLC improves CJ-speed error relative to LLF
and the wave propagates with clean HLLC/chemistry diagnostics. It is not final
detonation validation because CJ error remains above `5%`, post-front state
errors remain sizable, and HLLC planar dt/grid sensitivity has not been
completed.

### Current HLLC verification conclusion

**B. PARTIAL PASS** for the full HLLC 1D physics gates.

Allowed next work:

- refine HLLC reactive-shock front/reaction detection
- run controlled HLLC planar dt/grid sweep if needed

Still not allowed from this evidence alone:

- HLLC-P implementation
- AUSM+UP2 implementation
- 2D reactive detonation
- NS/diffusion/Soret/AMR/DLB
- claiming a fully validated detonation solver
