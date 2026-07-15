#!/bin/bash -l

#PBS -N pair_tempo_native
#PBS -A P19010000
#PBS -q casper
#PBS -l select=1:ncpus=1:mem=150GB
#PBS -l walltime=12:00:00
#PBS -J 1-30%6
#PBS -j oe
#PBS -o pair_tempo_native.log

# TEMPO-only native-resolution city pairing, one array element per June day.
# RUN (required) and CITY (optional) come from the qsub -v below.
# Submit with, e.g. (ATL + MEX on mxcat; ATL/DFW/LA/DEN on the CONUS runs):
#
#   for RUN in nonbiog biog grapes; do
#     for CITY in atl dfw la den; do
#       qsub -N pn_${RUN}_${CITY} -o pn_${RUN}_${CITY}.log \
#            -v RUN=$RUN,CITY=$CITY submit_pair_tempo_native.sh
#     done
#   done
#   for CITY in atl mex; do
#     qsub -N pn_mxcat_${CITY} -o pn_mxcat_${CITY}.log \
#          -v RUN=mxcat,CITY=$CITY submit_pair_tempo_native.sh
#   done
#
# Omit CITY to pair all six boxes in one (longer) job. Override resolution
# with -v RUN=...,CITY=...,OBS_GRID_RES=0.05 if 0.03 is too fine/slow.

cd /glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/hcho_ben_gaubert_cs/batch_scripts

module load conda
conda activate melodies-monet

export YMD=$(printf "202406%02d" ${PBS_ARRAY_INDEX})   # 20240601 .. 20240630

python pair_tempo_native.py