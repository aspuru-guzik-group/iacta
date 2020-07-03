#!/bin/bash
#SBATCH -A rrg-aspuru
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=8:00:00
#SBATCH --job-name minigabe-restart

module load NiaEnv/2018a
module load intel
module load python/3.6.4-anaconda5.1.0
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH=$HOME/xtb/bin:$PATH
export LOCALSCRATCH=$SLURM_TMPDIR
export MINIGABE=$HOME/ts-search

source activate ts-search
python3.7 $MINIGABE/rsearch.py -t 40 $@
