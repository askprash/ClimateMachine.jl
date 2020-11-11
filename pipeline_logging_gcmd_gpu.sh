#!/bin/bash

#SBATCH --nodes=1
#SBATCH --tasks-per-node=2 # number of MPI ranks per node
#SBATCH --gres=gpu:2       # GPUs per node; should equal tasks-per-node
#SBATCH --time=01:00:00

# Kill the job if anything fails
#set -euo pipefail
set -x # echo script

module purge;
#module load julia/1.4.2 hdf5/1.10.1 netcdf-c/4.6.1 cuda/10.0 openmpi/4.0.3_cuda-10.0 # CUDA-aware MPI
module load julia/1.5.2 hdf5/1.10.1 netcdf-c/4.6.1 cuda/10.0 openmpi/4.0.3_cuda-10.0 # CUDA-aware MPI

#export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:=1}
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:=1}
export JULIA_MPI_BINARY=system
export JULIA_CUDA_USE_BINARYBUILDER=false

# Import helper functions for this script
source ./helper_mod.sh

# User envirnoment setup
RUNNAME="gcm_spectrum"

# Change if CLIMA and VizCLIMA not saved in $HOME
CLIMA_HOME='/central/groups/esm/lenka/ClimateMachine.jl'

CLIMA_OUTPUT='/central/scratch/elencz/output/'$RUNNAME

# Choose CLIMA experiment script and VizCLIMA script
EXPERIMENT=$CLIMA_HOME'/experiments/AtmosGCM/GCMDriver/GCMDriver.jl'


directory_structure $CLIMA_OUTPUT
julia --project=$CLIMA_HOME -e 'using Pkg; Pkg.instantiate()'
julia --project=$CLIMA_HOME -e 'using Pkg; Pkg.precompile()'

write_into_runfile_from_list "$line" "$RUNNAME" "$EXPERIMENT" CLIMA_RUNFILE

# run a --project=$CLIMA_HOME -e ''each experiment listed in EXP_PARAM_FILE
# Run climate model
mpiexec julia --project=$CLIMA_HOME $CLIMA_RUNFILE --experiment=heldsuarez --diagnostics 0.1shours --monitor-courant-numbers 6shours --output-dir $CLIMA_NETCDF --checkpoint-at-end --checkpoint-dir $CLIMA_RESTART --init-moisture-profile zero --checkpoint 6shours --surface-flux bulk 
