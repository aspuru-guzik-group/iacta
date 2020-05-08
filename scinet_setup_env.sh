#!/bin/bash

module load intel
module load python/3.6.4-anaconda5.1.0
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH=$HOME/xtb/bin:$PATH
export LOCALSCRATCH=$SLURM_TMPDIR
export MINIGABE=$HOME/ts-search

source activate ts-search
