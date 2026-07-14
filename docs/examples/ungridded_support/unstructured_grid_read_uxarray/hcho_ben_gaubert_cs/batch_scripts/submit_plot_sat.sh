#!/bin/bash -l

#PBS -N plot_sat_data
#PBS -A P19010000
#PBS -q casper
#PBS -l select=1:ncpus=1:mem=48GB
#PBS -l walltime=04:00:00
#PBS -j oe
#PBS -o plot_sat.log

# One-time sat plotting job 

# Sat plotting job. Default 48GB fits one PLOT_GRID x PLOT_ONLY split;
# override per submission with qsub -l (see the loops below).

cd /glade/u/home/lcthompson/mm/MELODIES-MONET/docs/examples/ungridded_support/unstructured_grid_read_uxarray/hcho_ben_gaubert_cs/batch_scripts

module load conda
conda activate melodies-monet   

python plot_sat.py





# # standard runs obs-grid groups (timeseries/box/overlay on the 0.1deg pairs) 
# for RUN in nonbiog biog grapes mxcat; do
#   for ONLY in grp1_ grp1b_ grp4_ grp5; do
#     export MM_CONTROL=$PWD/control_${RUN}.yaml PLOT_GRID=obs PLOT_ONLY=$ONLY
#     qsub -N po_${RUN}_${ONLY%_} -o po_${RUN}_${ONLY%_}.log \
#          -l select=1:ncpus=1:mem=32GB -V submit_plot_sat.sh
#   done
# done

# # standard runs: model-grid groups (incl. diurnal/windowed/multibox) 
# for RUN in nonbiog biog grapes mxcat; do
#   for ONLY in grp1_ grp4_ grp5 grp6 grp7 grp8; do
#     export MM_CONTROL=$PWD/control_${RUN}.yaml PLOT_GRID=model PLOT_ONLY=$ONLY
#     qsub -N pm_${RUN}_${ONLY} -o pm_${RUN}_${ONLY}.log \
#          -l select=1:ncpus=1:mem=48GB -V submit_plot_sat.sh
#   done
# done

# # conservative model-space runs (TEMPO mesh products) 
# for RUN in nonbiog_cons biog_cons grapes_cons mxcat_cons; do
#   for ONLY in grp5 grp6 grp7 grp8; do
#     export MM_CONTROL=$PWD/control_${RUN}.yaml PLOT_GRID=model PLOT_ONLY=$ONLY
#     qsub -N pc_${RUN}_${ONLY} -o pc_${RUN}_${ONLY}.log \
#          -l select=1:ncpus=1:mem=48GB -V submit_plot_sat.sh
#   done
# done

# # native city products: obsgrid maps + swath footprints, both instruments 
# for RUN in nonbiog biog grapes mxcat; do
#   for CTL in control_${RUN}_native_cons.yaml control_${RUN}_tropomi_native_cons.yaml; do
#     [ -f "$CTL" ] || continue
#     export MM_CONTROL=$PWD/$CTL PLOT_GRID=obs PLOT_ONLY=grp_native
#     qsub -N pn_${CTL%.yaml} -o pn_${CTL%.yaml}.log \
#          -l select=1:ncpus=1:mem=24GB -V submit_plot_sat.sh
#   done
#   for CTL in control_${RUN}_native_cons_swath.yaml control_${RUN}_tropomi_native_cons_swath.yaml; do
#     [ -f "$CTL" ] || continue
#     export MM_CONTROL=$PWD/$CTL PLOT_GRID=obs PLOT_ONLY=grp_swath
#     qsub -N ps_${CTL%.yaml} -o ps_${CTL%.yaml}.log \
#          -l select=1:ncpus=1:mem=24GB -V submit_plot_sat.sh
#   done
# done




    # for GRID in obs model; do
    #   for ONLY in grp1_ grp1b_ grp2_ grp3 grp4_ grp5; do
    #     export PLOT_GRID=$GRID PLOT_ONLY=$ONLY
    #     tag=${GRID}_${ONLY%_}
    #     qsub -N p_${tag} -o p_${tag}.log -l select=1:ncpus=1:mem=48GB -V submit_plot_sat.sh
    #   done
    # done
    
# qsub -N p_co -o p_co.log -l select=1:ncpus=1:mem=48GB \
#      -v PLOT_ONLY=_co submit_plot_sat.sh