#!/bin/bash -l

#SBATCH -J test
#SBATCH -A m224
#SBATCH -p debug
#SBATCH -o test.out
#SBATCH -e test.err
#SBATCH -N 16
#SBATCH -t 0:25:00
##SBATCH --mail-type=ALL
#SBATCH -V
#SBATCH -C haswell

export I_MPI_FABRICS=ofi
export I_MPI_OFI_PROVIDER=gni
export I_MPI_OFI_LIBRARY=/global/common/cori/software/libfabric/1.6.1/gnu/lib/libfabric.so
export I_MPI_PMI_LIBRARY=/usr/lib64/slurmpmi/libpmi.so

cd $SLURM_SUBMIT_DIR
EXEPATH="./viriato"
ARGS=" 3em4_04 " 
NPROCS=$(expr $SLURM_JOB_NUM_NODES \* $SLURM_CPUS_ON_NODE / 2)

srun -n 512 $EXEPATH $ARGS
