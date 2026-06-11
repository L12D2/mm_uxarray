#!/bin/bash -l

#PBS -N pair_obs_data
#PBS -A P19010000
#PBS -q casper
#PBS -l select=1:ncpus=4:mem=120GB
#PBS -l walltime=02:00:00
#PBS -j oe
#PBS -o pair_dc8.log

# One-time DC8 aircraft pairing job 

cd /glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/asiaaq_cs_06082026/write_met_data/batch_scripts/

module load conda
conda activate melodies-monet   

python pair_aircraft.py