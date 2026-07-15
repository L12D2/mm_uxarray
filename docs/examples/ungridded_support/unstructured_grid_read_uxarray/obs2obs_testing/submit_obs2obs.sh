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


# for ONLY in em_no2col em_no2sfc em_o3sfc em_hchocol diurnal; do
#   export OBS2OBS_ONLY=$ONLY
#   qsub -N o2o_$ONLY -o o2o_$ONLY.log -l select=1:ncpus=1:mem=48GB -V submit_obs2obs.sh
# done

# for ONLY in em_no2col em_hchocol diurnal; do
#   export OBS2OBS_ONLY=$ONLY
#   qsub -N o2o_$ONLY -o o2o_$ONLY.log -l select=1:ncpus=1:mem=96GB -V submit_obs2obs.sh
# done

# for ONLY in em_no2col; do
#   export OBS2OBS_ONLY=$ONLY
#   qsub -N o2o_$ONLY -o o2o_$ONLY.log -l select=1:ncpus=1:mem=128GB -V submit_obs2obs.sh
# done

# # these ones require more mem 
# for ONLY in coupling operator grid no2_bias; do
#   export OBS2OBS_ONLY=$ONLY
#   qsub -N o2o_$ONLY -o o2o_$ONLY.log -l select=1:ncpus=1:mem=200GB -V submit_obs2obs.sh
# done



# for ONLY in no2_bias; do
#   export OBS2OBS_ONLY=$ONLY
#   qsub -N o2o_$ONLY -o o2o_$ONLY.log -l select=1:ncpus=1:mem=200GB -V submit_obs2obs.sh
# done

