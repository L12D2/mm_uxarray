#!/bin/bash -l

#PBS -N plot_sat_data
#PBS -A P19010000
#PBS -q casper
#PBS -l select=1:ncpus=1:mem=200GB
#PBS -l walltime=04:00:00
#PBS -j oe
#PBS -o plot_sat.log

# One-time sat plotting job 

cd /glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/hcho_ben_gaubert_cs/batch_scripts

module load conda
conda activate melodies-monet   

python plot_sat-mxcat.py

    # for GRID in obs model; do
    #   for ONLY in grp1_ grp1b_ grp2_ grp3 grp4_ grp5; do
    #     export PLOT_GRID=$GRID PLOT_ONLY=$ONLY
    #     tag=${GRID}_${ONLY%_}
    #     qsub -N p_${tag} -o p_${tag}.log -l select=1:ncpus=1:mem=48GB -V submit_plot_sat-bio.sh
    #   done
    # done
