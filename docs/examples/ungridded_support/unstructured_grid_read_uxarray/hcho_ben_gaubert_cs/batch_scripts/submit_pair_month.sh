#!/bin/bash -l

#PBS -N pair_month_sat_data
#PBS -A P19010000
#PBS -q casper
#PBS -l select=1:ncpus=1:mem=150GB
#PBS -l walltime=06:00:00
#PBS -J 1-30%6
#PBS -j oe
#PBS -o pair_month_sat_data.log

# One month of standard CONUS sat pairing (radius_mean products), one array
# element per June day. RUN selects the emissions run (see RUNS in
# pair_daily_sat.py); OBS_GROUP optionally splits by obs product.
# Peak memory ~74 GB per element (full model day) -> 96 GB.
#
#   for RUN in nonbiog biog grapes mxcat; do
#     qsub -N ptm_${RUN} -o ptm_${RUN}.log -v RUN=$RUN submit_pair_month.sh
#   done
#
#   # or split TROPOMI products into their own jobs:
#   # qsub -N ptm_${RUN}_no2 -v RUN=$RUN,OBS_GROUP=tropomi_l2_no2 submit_pair_month.sh


cd /glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/hcho_ben_gaubert_cs/batch_scripts

module load conda
conda activate melodies-monet   

export YMD=$(printf "202406%02d" ${PBS_ARRAY_INDEX})   # 20240601 .. 20240630
#export YMD=$(printf "202403%02d" ${PBS_ARRAY_INDEX})   # 20240301 .. 20240331

python pair_daily_sat.py

# for P in tropomi_l2_no2 tropomi_l2_hcho tropomi_l2_co; do
#   qsub -N pair_${P} -o pair_${P}.log \
#        -l select=1:ncpus=1:mem=96GB \
#        -v OBS_GROUP=$P \
#        submit_pair_month.sh 
# done


# for RUN in nonbiog biog grapes mxcat; do
#   for P in tropomi_l2_no2 tropomi_l2_hcho tropomi_l2_co; do
#     qsub -N ptm_${RUN}_${P} -o ptm_${RUN}_${P}.log \
#          -v RUN=$RUN,OBS_GROUP=$P \
#          submit_pair_month.sh
#   done
# done

# for RUN in nonbiog biog grapes mxcat; do
#   for P in tempo_l2_no2 tempo_l2_hcho; do
#     qsub -N ptm_${RUN}_${P} -o ptm_${RUN}_${P}.log \
#          -v RUN=$RUN,OBS_GROUP=$P \
#          submit_pair_month.sh
#   done
# done