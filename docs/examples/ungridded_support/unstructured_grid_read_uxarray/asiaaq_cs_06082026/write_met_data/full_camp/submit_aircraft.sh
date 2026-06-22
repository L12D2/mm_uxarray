#!/bin/bash -l

#PBS -N pair_dc8_data
#PBS -A P19010000
#PBS -q casper
#PBS -l select=1:ncpus=4:mem=400GB
#PBS -l walltime=10:00:00
#PBS -o pair_dc8_full.log

# One-time DC8 aircraft pairing job 

cd /glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/asiaaq_cs_06082026/write_met_data/full_camp/

module load conda
conda activate melodies-monet   

python campaign_driver.py 