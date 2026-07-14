#!/bin/bash -l

#PBS -N obs2obstesting
#PBS -A P19010000
#PBS -q casper
#PBS -l select=1:ncpus=1:mem=96GB
#PBS -l walltime=02:00:00
#PBS -j oe
#PBS -o obs2obs_testing.log

# One-time sat pairing job 

cd /glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/obs2obs_testing

module load conda
conda activate melodies-monet   

python obs2obs.py

# export PLOT_GRID=obs
# qsub -N plot_obs   -o plot_obs.log   -l select=1:ncpus=1:mem=96GB -V submit_sat.sh
# export PLOT_GRID=model
# qsub -N plot_model -o plot_model.log -l select=1:ncpus=1:mem=96GB -V submit_sat.sh

# for ONLY in em_no2col em_no2sfc em_o3sfc em_hchocol no2_bias coupling skill diurnal multiscale; do
#   export OBS2OBS_ONLY=$ONLY
#   qsub -N o2o_$ONLY -o o2o_$ONLY.log -l select=1:ncpus=1:mem=48GB -V submit_obs2obs.sh
# done