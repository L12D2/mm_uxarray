#!/bin/bash -l

#PBS -N plot_test
#PBS -A P19010000
#PBS -q casper
#PBS -l select=1:ncpus=1:mem=50GB
#PBS -l walltime=02:00:00
#PBS -j oe
#PBS -o plot_test.log

# One-time sat plotting job 

cd /glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/hcho_ben_gaubert_cs/batch_scripts

module load conda
conda activate melodies-monet   

python plot_test.py
