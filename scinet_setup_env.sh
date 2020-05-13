#!/bin/bash

# This nasty piece of bash script yields the current directory where this file resides..
pushd . > /dev/null
SCRIPT_PATH="${BASH_SOURCE[0]}"
if ([ -h "${SCRIPT_PATH}" ]); then
	while([ -h "${SCRIPT_PATH}" ]); do cd `dirname "$SCRIPT_PATH"`; 
		SCRIPT_PATH=`readlink "${SCRIPT_PATH}"`; done
fi
cd `dirname ${SCRIPT_PATH}` > /dev/null
SCRIPT_PATH=`pwd`;
popd  > /dev/null

module load intel
module load python/3.6.4-anaconda5.1.0
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH=$HOME/xtb/bin:$PATH
export LOCALSCRATCH=$SLURM_TMPDIR
export MINIGABE="$SCRIPT_PATH"
echo "Minigabe in "$MINIGABE""

source activate ts-search
