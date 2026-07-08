#!/bin/bash -l
#PBS -N repair_time
#PBS -A P19010000
#PBS -q casper
#PBS -l select=1:ncpus=1:mem=96GB
#PBS -l walltime=06:00:00
#PBS -j oe
#PBS -o repair_time.log
cd /glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/hcho_ben_gaubert_cs/batch_scripts
module load conda
conda activate melodies-monet
line=$(sed -n "${PBS_ARRAY_INDEX}p" repair_time_manifest.txt)
[ -z "$line" ] && exit 0
read -r CFG PROD YMD <<< "$line"
case "$CFG" in
  nonbiog_refera5_dust) SCRIPT=pair_daily_sat.py ;;
  biog_refera5_dust)    SCRIPT=bio_pair_daily_sat.py ;;
  mxcat)                SCRIPT=mxcat_pair_daily_sat.py ;;
  *) echo "unknown config $CFG"; exit 1 ;;
esac
export YMD OBS_GROUP="$PROD"
echo "repair[${PBS_ARRAY_INDEX}] $CFG $PROD $YMD"
python "$SCRIPT"
