#!/bin/bash -l

#PBS -N pair_tropomi_native
#PBS -A P19010000
#PBS -q casper
#PBS -l select=1:ncpus=1:mem=150GB
#PBS -l walltime=12:00:00
#PBS -J 1-30%6
#PBS -j oe
#PBS -o pair_tropomi_native.log

# Native-resolution TROPOMI city pairing, one array element per June day.
# RUN (required) and CITY/REGRID_TARGET/TROPOMI_PRODUCTS (optional) come from
# the qsub -v below. Submit with, e.g.:
#
#   for RUN in nonbiog biog grapes; do
#     for CITY in atl dfw la den; do
#       qsub -N tn_${RUN}_${CITY} -o tn_${RUN}_${CITY}.log \
#            -v RUN=$RUN,CITY=$CITY,REGRID_METHOD=conservative,REGRID_TARGET=swath \
#            submit_pair_tropomi_native.sh
#     done
#   done
#   for CITY in atl mex; do
#     qsub -N tn_mxcat_${CITY} -o tn_mxcat_${CITY}.log \
#          -v RUN=mxcat,CITY=$CITY,REGRID_METHOD=conservative,REGRID_TARGET=swath \
#          submit_pair_tropomi_native.sh
#   done
#
# Restrict products with -v ...,TROPOMI_PRODUCTS=no2,hcho  (default: no2,hcho,co).

cd /glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/hcho_ben_gaubert_cs/batch_scripts

module load conda
conda activate melodies-monet

export YMD=$(printf "202406%02d" ${PBS_ARRAY_INDEX})   # 20240601 .. 20240630

python pair_tropomi_native.py