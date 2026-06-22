#!/bin/bash -l

#PBS -N pair_sfc_data
#PBS -A P19010000
#PBS -q casper
#PBS -l select=1:ncpus=4:mem=300GB
#PBS -l walltime=06:00:00
#PBS -j oe
#PBS -o pair_sfc_full.log

# One-time sfc pairing job 

cd /glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/asiaaq_cs_06082026/write_met_data/full_camp/

module load conda
conda activate melodies-monet   

python sfc_driver.py
