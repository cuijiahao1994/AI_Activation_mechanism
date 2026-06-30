#!/bin/bash
#SBATCH -N 1
#SBATCH -n 56
#SBATCH --ntasks-per-node=56
#SBATCH --partition=E116
#SBATCH --time=02:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err

source /srv/nfs/share/opt/intel/oneapi-2025.0/setvars.sh
export PATH=/srv/nfs/share/app/vasp.6.4.3/bin:$PATH
ulimit -s unlimited
mpirun -np $SLURM_NPROCS vasp_gam
