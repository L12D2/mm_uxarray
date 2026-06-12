
#!/bin/bash -l
#PBS -N tropomi_co_cesmse
#PBS -A P19010000
#PBS -q casper
#PBS -l select=1:ncpus=8:mem=200GB
#PBS -l walltime=04:00:00
#PBS -j oe
#PBS -o tropomi_co_cesmse.log

cd /glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/tropomi_testing/
mkdir -p output_co

module load conda
conda activate /glade/work/lcthompson/conda-envs/melodies-monet

python run_tropomi_co.py