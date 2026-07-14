#!/bin/bash -l

#PBS -N plot_sfc_data
#PBS -A P19010000
#PBS -q casper
#PBS -l select=1:ncpus=1:mem=80GB
#PBS -l walltime=04:00:00
#PBS -j oe
#PBS -o plot_sfc.log

# One-time sat plotting job 

  # for RUN in nonbiog biog grapes mxcat; do
  #   qsub -N psfc_${RUN} -o psfc_${RUN}.log -v RUN=$RUN submit_plot_sfc.sh
  # done

cd /glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/hcho_ben_gaubert_cs/batch_scripts

module load conda
conda activate melodies-monet   

python pair_daily_sfc.py
