# SU2_deto Server Build Report

## Purpose
This report records the reproducible server-side configuration, build, and validation workflow for:

- `SU2_deto`: `/home/jmyang/detonationFoam/SU2_deto`
- `Cantera`: `/home/jmyang/detonationFoam/third_party/install/cantera-3.2.0`
- `SUNDIALS`: `/home/jmyang/detonationFoam/third_party/install/sundials-7.7.0`

The goal was to enable, in one consistent server environment:

- Intel MPI parallel build/runtime
- Cantera
- SUNDIALS/CVODE
- Mutation++ / MPP
- `mpirun -np 2 SU2_CFD ...`

## Key Diagnosis Before Any Source Change

### 1. MPI runtime and wrapper status
The server already had Intel MPI available and working:

- `mpirun`: `/opt/intel/oneapi/mpi/2021.16/bin/mpirun`
- `mpicc`: `/opt/intel/oneapi/mpi/2021.16/bin/mpicc`
- `mpicxx`: `/opt/intel/oneapi/mpi/2021.16/bin/mpicxx`

A standalone MPI hello-world test compiled and ran successfully with 2 ranks after running outside the local sandbox restriction.

### 2. Cantera discovery
Meson and `pkg-config` both found Cantera correctly:

- `Run-time dependency cantera found: YES 3.2.0`
- `cantera.pc` exists under `.../lib/pkgconfig`

### 3. SUNDIALS discovery
Meson did **not** report the individual SUNDIALS dependencies as `YES` through `pkg-config` / direct CMake dependency discovery:

- `sundials_cvode found: NO`
- `sundials_nvecserial found: NO`
- `sundials_core found: NO`
- `sundials_sunmatrixdense found: NO`
- `sundials_sunlinsoldense found: NO`

However, this is **not** a blocker for this SU2 tree, because `SU2_deto/meson.build` already contains a fallback path controlled by:

- `-Denable-sundials=true`
- `-Dsundials_root=<path>`

That fallback manually injects the SUNDIALS include path and link flags from `SUNDIALS_ROOT`.

### 4. Why MPI dependency detection still failed even with Intel MPI wrappers
The old Meson log showed that `dependency('mpi')` was still failing even when Intel wrapper compilers had been tried.

Instead of patching wrappers again, the existing built-in SU2 option was used:

- `-Dcustom-mpi=true`

This is already supported by `SU2_deto/meson.build` and tells Meson to trust the explicitly selected MPI wrapper compilers (`CC=mpicc`, `CXX=mpicxx`) rather than requiring `dependency('mpi')` to succeed.

This is the smallest clean solution for this server.

## No Source Patch Was Required
No tracked SU2 source file needed to be patched for this server build.

Important note: the copied server directory is **not** a valid full git checkout. `/home/jmyang/detonationFoam/.git` exists but does not contain a usable repository state, so `git diff` is not available on the server copy.

Therefore:

- no source patch was applied in this session;
- no tracked-file `git diff` can be generated from this packaged server copy.

## Reproducible Environment Script
The server environment script created for this build is:

- `/home/jmyang/detonationFoam/SU2_deto/env_su2_deto.sh`

It sets:

- `MPI_HOME=/opt/intel/oneapi/mpi/2021.16`
- Intel MPI `bin` first on `PATH`
- the Conda environment containing `mpi4py`
- `CANTERA_ROOT`
- `SUNDIALS_ROOT`
- `LD_LIBRARY_PATH`
- `PKG_CONFIG_PATH`
- `CMAKE_PREFIX_PATH`
- `MPP_DATA_DIRECTORY`
- `CC` / `CXX` to Intel MPI wrappers

## Final Configure Command
The successful clean configuration command was:

```bash
source /home/jmyang/detonationFoam/SU2_deto/env_su2_deto.sh
cd /home/jmyang/detonationFoam/SU2_deto
rm -rf build install
./meson.py setup build \
  -Denable-pywrapper=true \
  -Denable-mpp=true \
  -Denable-sundials=true \
  -Dsundials_root=$SUNDIALS_ROOT \
  -Denable-cantera=true \
  -Dcantera_root=$CANTERA_ROOT \
  -Dcustom-mpi=true \
  --prefix=$PWD/install
```

## Why `custom-mpi=true` Was Necessary
Meson 0.61 on this server did not reliably resolve Intel MPI through `dependency('mpi')`, even though:

- Intel MPI wrappers worked;
- MPI hello-world compiled;
- MPI hello-world ran correctly.

`custom-mpi=true` is therefore a **server adaptation using an existing SU2 build option**, not a local hack. It avoids mixing OpenMPI and Intel MPI, and it keeps the build on:

- `CC=/opt/intel/oneapi/mpi/2021.16/bin/mpicc`
- `CXX=/opt/intel/oneapi/mpi/2021.16/bin/mpicxx`

No `meson.build` patch was required.

## Build Result
The build completed successfully with:

```bash
./ninja -C build install
```

Installed binary path:

- `/home/jmyang/detonationFoam/SU2_deto/install/bin/SU2_CFD`

Other installed components include:

- `SU2_DEF`
- `SU2_DOT`
- `SU2_GEO`
- `SU2_SOL`
- `_pysu2.so`

## Meson Audit Evidence
The copied Meson log is:

- `/home/jmyang/detonationFoam/SU2_deto/build_audit/meson-log.txt`

Key lines:

- `C compiler for the host machine: /opt/intel/oneapi/mpi/2021.16/bin/mpicc`
- `C++ compiler for the host machine: /opt/intel/oneapi/mpi/2021.16/bin/mpicxx`
- `Run-time dependency cantera found: YES 3.2.0`
- `Run-time dependency sundials_* found: NO (tried pkgconfig and cmake)`
- `Message: Using mpi4py from /home/jmyang/miniconda3/envs/su2_deto/lib/python3.10/site-packages/mpi4py/include`
- `Subproject Mutationpp finished.`

Interpretation:

- Cantera was discovered directly.
- SUNDIALS was linked through the `sundials_root` fallback path.
- MPI was enabled through wrapper compilers plus `custom-mpi=true`.
- Mutation++ / MPP was enabled and built.

## Installed Binary Validation
After sourcing `env_su2_deto.sh`, the following checks succeeded:

```bash
which SU2_CFD
SU2_CFD -h | head -50
mpirun -np 2 SU2_CFD -h | head -50
```

### Relevant `ldd` lines
`ldd $(which SU2_CFD)` showed:

- `libmutation__.so => /home/jmyang/detonationFoam/SU2_deto/build/subprojects/Mutationpp/libmutation__.so`
- `libcantera.so.3 => /home/jmyang/detonationFoam/third_party/install/cantera-3.2.0/lib/libcantera.so.3`
- `libsundials_cvode.so.7 => /home/jmyang/detonationFoam/third_party/install/sundials-7.7.0/lib/libsundials_cvode.so.7`
- `libsundials_core.so.7 => /home/jmyang/detonationFoam/third_party/install/sundials-7.7.0/lib/libsundials_core.so.7`
- `libmpi.so.12 => /opt/intel/oneapi/mpi/2021.16/lib/libmpi.so.12`
- `libsundials_cvodes.so.7 => /home/jmyang/detonationFoam/third_party/install/sundials-7.7.0/lib/libsundials_cvodes.so.7`
- `libsundials_idas.so.6 => /home/jmyang/detonationFoam/third_party/install/sundials-7.7.0/lib/libsundials_idas.so.6`

No `not found` entries appeared in the filtered output.

## Minimal Runtime Verification
Two small cases were run under Intel MPI with 2 ranks.

### 1. Minimal inert case
Command:

```bash
source /home/jmyang/detonationFoam/SU2_deto/env_su2_deto.sh
cd /home/jmyang/detonationFoam/SU2_test_case/verified_detonation_euler_cases/cases/03_inert_uniform_freestream
mpirun -np 2 SU2_CFD uniform_freestream.cfg
```

Result:

- Success: yes
- MPI ranks: 2
- Chemistry/CVODE error: no
- Output files observed:
  - `history_uniform.csv`
  - `restart_uniform_00019.csv`

### 2. Minimal DETONATION_EULER chemistry-only case
Command:

```bash
source /home/jmyang/detonationFoam/SU2_deto/env_su2_deto.sh
cd /home/jmyang/detonationFoam/SU2_test_case/verified_detonation_euler_cases/cases/00_chemistry_only_0d
mpirun -np 2 SU2_CFD chemistry_only.cfg
```

Result:

- Success: yes
- MPI ranks: 2
- Cantera mechanism load: yes
- Chemistry/CVODE fatal error: no
- Output files observed:
  - `history_chemistry_only.csv`
  - `restart_chemistry_only_00049.csv`

## Final Enablement Status
The current server build has truly enabled:

- MPI: **YES**
- Cantera: **YES**
- SUNDIALS / CVODE: **YES**
- Mutation++ / MPP: **YES**

## How to Run with More Ranks
The same runtime path can be used with more ranks, for example:

```bash
source /home/jmyang/detonationFoam/SU2_deto/env_su2_deto.sh
cd /path/to/case
mpirun -np 8 SU2_CFD your_case.cfg
```

For a scheduler-managed server, use the scheduler to allocate ranks and then launch the same executable under the allocated MPI environment. The critical requirement is to keep using the same Intel MPI environment from `env_su2_deto.sh` and not fall back to `/usr/bin/mpirun` or `/usr/bin/mpic++` from OpenMPI.

## Remaining Server Caveat
The top-level copied project directory is not a normal git checkout on the server, so future server-side source patches should be tracked carefully in a local report file unless the repository is recopied with full git metadata.

## Slurm Submission Validation

A generic Slurm submission wrapper was validated on the server:

- wrapper script: `/home/jmyang/detonationFoam/sub_deto.sh`
- environment wrapper: `/home/jmyang/detonationFoam/env_su2_deto.sh`

### Minimal fixes required for Slurm
Two small server-side fixes were needed to make `sbatch` reliable:

1. The top-level `env_su2_deto.sh` was converted into a compatibility wrapper that simply sources:
   - `/home/jmyang/detonationFoam/SU2_deto/env_su2_deto.sh`

2. `SU2_deto/env_su2_deto.sh` was made safe under `set -u`.
   The root cause of the first failed Slurm job was that Intel oneAPI `setvars.sh` references `OCL_ICD_FILENAMES` as an unset variable in batch shells running with `nounset` enabled. The fix temporarily disables `nounset` only around:
   - `source /opt/intel/oneapi/setvars.sh`

### Generic Slurm usage
The validated usage is:

```bash
cd /path/to/case_directory
sbatch /home/jmyang/detonationFoam/sub_deto.sh your_case.cfg
```

The wrapper:

- changes into the case directory,
- sources the server SU2 environment,
- launches:
  - `srun -n $SLURM_NTASKS SU2_CFD <cfg>`

### Slurm validation jobs completed
Three small jobs were submitted **serially** and monitored to completion.

1. Inert uniform freestream
   - Job ID: `2560`
   - State: `COMPLETED`
   - Command:
     ```bash
     cd /home/jmyang/detonationFoam/SU2_test_case/verified_detonation_euler_cases/cases/03_inert_uniform_freestream
     sbatch /home/jmyang/detonationFoam/sub_deto.sh uniform_freestream.cfg
     ```

2. DETONATION_EULER chemistry-only
   - Job ID: `2561`
   - State: `COMPLETED`
   - Command:
     ```bash
     cd /home/jmyang/detonationFoam/SU2_test_case/verified_detonation_euler_cases/cases/00_chemistry_only_0d
     sbatch /home/jmyang/detonationFoam/sub_deto.sh chemistry_only.cfg
     ```

3. DETONATION_EULER closed flow + chemistry
   - Job ID: `2562`
   - State: `COMPLETED`
   - Command:
     ```bash
     cd /home/jmyang/detonationFoam/SU2_test_case/verified_detonation_euler_cases/cases/01_closed_flow_plus_chemistry_0d
     sbatch /home/jmyang/detonationFoam/sub_deto.sh flow_plus_chemistry.cfg
     ```

### Practical note
For larger production-style runs, keep the same wrapper and adjust Slurm resources at submit time, for example:

```bash
sbatch --ntasks-per-node=8 /home/jmyang/detonationFoam/sub_deto.sh your_case.cfg
```

This preserves the same Intel MPI + Cantera + SUNDIALS + Mutation++ environment while scaling the rank count.

## Standard Daily Workflow (Simple Version)

### One command you should always run first in a new terminal
For normal daily use, the standard entry point is:

```bash
source /home/jmyang/detonationFoam/env_su2_deto.sh
```

This is the recommended **single** setup step after opening a new terminal.

You do **not** need to separately run both of the following every time:

- `module load mpi/2021.16`
- `conda activate su2_deto`

The environment wrapper already prepares the required MPI, Python, Cantera, SUNDIALS, and Mutation++ runtime/build paths.

### What the environment wrapper provides
After sourcing `env_su2_deto.sh`, the session is prepared with:

- Intel MPI runtime and wrapper compiler paths
- the `su2_deto` Python environment on `PATH`
- `SU2_HOME`
- `SU2_RUN`
- `CANTERA_ROOT`
- `SUNDIALS_ROOT`
- `LD_LIBRARY_PATH`
- `PKG_CONFIG_PATH`
- `CMAKE_PREFIX_PATH`
- `MPP_DATA_DIRECTORY`
- `CC` / `CXX` set to Intel MPI wrappers

## Standard Development / Run Procedures

### A. Open a new terminal and prepare the environment
```bash
cd /home/jmyang/detonationFoam
source /home/jmyang/detonationFoam/env_su2_deto.sh
```

### B. Run a case directly with `mpirun`
```bash
source /home/jmyang/detonationFoam/env_su2_deto.sh
cd /path/to/case_directory
mpirun -np 2 SU2_CFD your_case.cfg
```

To use more ranks:
```bash
mpirun -np 8 SU2_CFD your_case.cfg
```

### C. Run a case through Slurm
```bash
cd /path/to/case_directory
sbatch /home/jmyang/detonationFoam/sub_deto.sh your_case.cfg
```

Or request more ranks at submit time:
```bash
sbatch --ntasks-per-node=8 /home/jmyang/detonationFoam/sub_deto.sh your_case.cfg
```

### D. Rebuild SU2_deto
```bash
source /home/jmyang/detonationFoam/env_su2_deto.sh
cd /home/jmyang/detonationFoam/SU2_deto

rm -rf build install
./meson.py setup build \
  -Denable-pywrapper=true \
  -Denable-mpp=true \
  -Denable-sundials=true \
  -Dsundials_root=$SUNDIALS_ROOT \
  -Denable-cantera=true \
  -Dcantera_root=$CANTERA_ROOT \
  -Dcustom-mpi=true \
  --prefix=$PWD/install

./ninja -C build install
```

### E. Quick verification after rebuild
```bash
source /home/jmyang/detonationFoam/env_su2_deto.sh
which SU2_CFD
SU2_CFD -h | head -20
mpirun -np 2 SU2_CFD -h | head -20
```

## Git / GitHub Workflow on the Server

### Server-side git baseline
A usable git repository now exists inside:

- `/home/jmyang/detonationFoam/SU2_deto`

This means future source changes can be inspected with:

```bash
cd /home/jmyang/detonationFoam/SU2_deto
git status
git diff
git log --oneline -5
```

### Recommended workflow before and after edits
```bash
cd /home/jmyang/detonationFoam/SU2_deto
git status
git diff
# make changes
git status
git diff
```

### GitHub push target
Target repository provided by the user:

- `git@github.com:jillywillitzer1966-source/SU2_75aaa.git`

### GitHub access test
A read-only access test succeeded with:

```bash
GIT_SSH_COMMAND='ssh -o BatchMode=yes -o StrictHostKeyChecking=accept-new -o ConnectTimeout=8' \
  git ls-remote git@github.com:jillywillitzer1966-source/SU2_75aaa.git HEAD
```

The command completed without authentication failure, which indicates that SSH access from this server to the target repository is available.
