#!/bin/bash -l

#PBS -N test_curtain
#PBS -A P19010000
#PBS -q casper
#PBS -l select=1:ncpus=4:mem=50GB
#PBS -l walltime=01:00:00
#PBS -j oe
#PBS -o test_curtain.log

# One-time DC8 aircraft pairing job 

cd /glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/readthedocs/

module load conda
conda activate melodies-monet   

python pair_aircraft.py