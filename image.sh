#!/bin/bash
#SBATCH -J dust_therm
#SBATCH -t 24:00:00
#SBATCH -N 1 -n 24

export OMP_NUM_THREADS=24
source ~/radmc_setup.sh

#radmc3d mctherm setthreads 24

for i in 0.1 1.0 10.0
do
  radmc3d image lambda $i nostar nphot_scat 10000000 setthreads 24 npix 256 incl 180 phi 0 sizeau 1
  mv image.out image${i}_nostar.out
  radmc3d image lambda $i inclstar nphot_scat 10000000 setthreads 24 npix 256 incl 180 phi 0 sizeau 1
  mv image.out image${i}_inclstar.out
done
