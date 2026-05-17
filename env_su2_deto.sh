#!/usr/bin/env bash
# Reproducible build/run environment for SU2_deto on the server.
# Usage:
#   source /home/jmyang/detonationFoam/SU2_deto/env_su2_deto.sh
#   cd /home/jmyang/detonationFoam/SU2_deto
#   ./meson.py setup build ...
#   ./ninja -C build install
#   mpirun -np 2 SU2_CFD -h

export DETO_ROOT=/home/jmyang/detonationFoam
export SU2_HOME=$DETO_ROOT/SU2_deto
export SU2_RUN=$SU2_HOME/install/bin
export SU2_TEST_CASE=$DETO_ROOT/SU2_test_case

export CANTERA_ROOT=$DETO_ROOT/third_party/install/cantera-3.2.0
export SUNDIALS_ROOT=$DETO_ROOT/third_party/install/sundials-7.7.0
export MPI_HOME=/opt/intel/oneapi/mpi/2021.16
export CONDA_SU2_ENV=/home/jmyang/miniconda3/envs/su2_deto

# oneAPI setvars can fail inside shells running with `set -u` because some
# vendor scripts reference unset variables. Preserve the caller's nounset mode
# while still allowing the vendor script to populate optional runtime entries.
__su2_restore_nounset=0
case $- in
  *u*) __su2_restore_nounset=1; set +u ;;
esac
if [ -f /opt/intel/oneapi/setvars.sh ]; then
  source /opt/intel/oneapi/setvars.sh >/tmp/su2_deto_setvars.log 2>&1 || true
fi
if [ "$__su2_restore_nounset" -eq 1 ]; then
  set -u
fi
unset __su2_restore_nounset

export PATH=$CONDA_SU2_ENV/bin:$MPI_HOME/bin:$SU2_RUN:$PATH
export PYTHONPATH=$SU2_RUN:${PYTHONPATH:-}

export PKG_CONFIG_PATH=$CANTERA_ROOT/lib/pkgconfig:$CANTERA_ROOT/lib64/pkgconfig:$SUNDIALS_ROOT/lib/pkgconfig:$SUNDIALS_ROOT/lib64/pkgconfig:${PKG_CONFIG_PATH:-}
export CMAKE_PREFIX_PATH=$CANTERA_ROOT:$SUNDIALS_ROOT:${CMAKE_PREFIX_PATH:-}
export LD_LIBRARY_PATH=$MPI_HOME/lib:$CANTERA_ROOT/lib:$CANTERA_ROOT/lib64:$SUNDIALS_ROOT/lib:$SUNDIALS_ROOT/lib64:$SU2_HOME/build/subprojects/Mutationpp:${LD_LIBRARY_PATH:-}

export MPP_DATA_DIRECTORY=$SU2_HOME/subprojects/Mutationpp/data
export MPLCONFIGDIR=$DETO_ROOT/mplconfig

# Build-time compiler selection. This keeps Meson on Intel MPI wrappers
# and avoids mixing OpenMPI wrappers from /usr/bin.
export CC=$MPI_HOME/bin/mpicc
export CXX=$MPI_HOME/bin/mpicxx
