#!/bin/bash
#SBATCH -J dust_therm
#SBATCH -t 24:00:00
#SBATCH -N 1 -n 24

export OMP_NUM_THREADS=24
source ~/radmc_setup.sh

radmc3d mctherm setthreads 24
