#!/bin/bash

#SBATCH -t 04:00:00
#SBATCH -n 12
#SBATCH --mem-per-cpu=1G

module load gcc openmpi/4.1.1_ft3 python 

export CPPFLAGS="  -O3"
export CFLAGS=" -O3"

export HDF5_BASE=/opt/cesga/2020/software/Compiler/gcccore/system/hdf5/1.12.1
export PYTHONPATH=/home/csic/gim/dro/.local/lib/python3.9/site-packages

echo "srun ./$1 $2 $3" 
srun ./$1 $2 $3
