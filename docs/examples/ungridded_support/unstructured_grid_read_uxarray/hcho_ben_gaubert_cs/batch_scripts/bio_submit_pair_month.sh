#!/bin/bash -l

#PBS -N pair_month_sat_data
#PBS -A P19010000
#PBS -q casper
#PBS -l select=1:ncpus=1:mem=180GB
#PBS -l walltime=02:00:00
#PBS -J 1-30%6
#PBS -j oe
#PBS -o pair_month_sat_data.log

# One month sat pairing 

cd /glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/hcho_ben_gaubert_cs/batch_scripts

module load conda
conda activate melodies-monet   

export YMD=$(printf "202406%02d" ${PBS_ARRAY_INDEX})   # 20240601 .. 20240630
#export YMD=$(printf "202403%02d" ${PBS_ARRAY_INDEX})   # 20240301 .. 20240331

python bio_pair_daily_sat.py

# for P in tempo_l2_hcho tempo_l2_no2 tropomi_l2_no2 tropomi_l2_hcho tropomi_l2_co; do
#   qsub -N pair_${P} -o pair_${P}.log \
#        -l select=1:ncpus=1:mem=96GB \
#        -v OBS_GROUP=$P \
#        bio_submit_pair_month.sh 
# done

# for P in tropomi_l2_no2 tropomi_l2_hcho tropomi_l2_co; do
#   qsub -N pair_${P} -o pair_${P}.log \
#        -l select=1:ncpus=1:mem=96GB \
#        -v OBS_GROUP=$P \
#        bio_submit_pair_month.sh 
# done