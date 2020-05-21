#!/bin/bash
#SBATCH -A rrg-aspuru
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=8:00:00
#SBATCH --job-name minigabe-restart

module load gnu-parallel 
module load intel
module load python/3.6.4-anaconda5.1.0
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH=$HOME/xtb/bin:$PATH
export LOCALSCRATCH=$SLURM_TMPDIR
export MINIGABE=$HOME/ts-search

source activate ts-search
cat $MINIGABE/test-set/all | parallel --joblog $(date +%d_%m_%Y).log --resume -j 1 "python $MINIGABE/rsearch-restart.py $MINIGABE/test-set/{}/user.yaml -o {} -w -t 40"
