#!/bin/bash -l

#PBS -N plot_air_data
#PBS -A P19010000
#PBS -q casper
#PBS -l select=1:ncpus=4:mem=200GB
#PBS -l walltime=06:00:00
#PBS -j oe
#PBS -o air_plot.log

# One-time air plot job 

cd /glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/asiaaq_cs_06082026/

module load conda
conda activate melodies-monet   

python air_plot_driver.py
