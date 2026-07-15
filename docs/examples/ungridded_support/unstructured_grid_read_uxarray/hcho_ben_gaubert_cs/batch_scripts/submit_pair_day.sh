#!/bin/bash -l

#PBS -N pair_day_sat
#PBS -A P19010000
#PBS -q casper
#PBS -l select=1:ncpus=1:mem=150GB
#PBS -l walltime=06:00:00
#PBS -j oe

# Re-pair a SINGLE day -- fills gaps left by failed/missing month-array
# elements. Unlike submit_pair_month.sh (YMD from PBS_ARRAY_INDEX), this takes
# YMD in via `qsub -v` so you can target only the missing dates.
# RUN selects the emissions run; OBS_GROUP optionally restricts obs products.
#
#   for RUN in grapes biog; do for DD in 07 08 14; do          # <- missing days
#     qsub -N rp_${RUN}_${DD} -o rp_${RUN}_${DD}.log \
#          -v RUN=$RUN,YMD=202406${DD},OBS_GROUP=tempo_l2_no2 submit_pair_day.sh
#   done; done

cd /glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/hcho_ben_gaubert_cs/batch_scripts

module load conda
conda activate melodies-monet

: "${YMD:?set YMD=202406DD via qsub -v}"      # fail loudly if unset (avoids 20240600)

python pair_daily_sat.py

#### ^ one day paiur