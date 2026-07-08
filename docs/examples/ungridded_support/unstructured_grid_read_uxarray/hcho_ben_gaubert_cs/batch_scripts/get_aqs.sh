#!/bin/bash -l

#PBS -N get_aqs
#PBS -A P19010000
#PBS -q casper
#PBS -l select=1:ncpus=1:mem=64GB
#PBS -l walltime=06:00:00
#PBS -j oe
#PBS -o get_aqs.log

# Download AQS surface obs for June 2024 (retry on failure)

cd /glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/hcho_ben_gaubert_cs/output

module load conda
conda activate melodies-monet

n=0
until melodies-monet get-aqs -s 2024-06-01 -e 2024-07-01 --no-compress; do
  n=$((n+1))
  if [ "$n" -ge 20 ]; then
    echo "get-aqs failed $n times; giving up." >&2
    exit 1
  fi
  echo "Download failed (attempt $n), retrying in 30s..." >&2
  sleep 30
done
echo "get-aqs done."