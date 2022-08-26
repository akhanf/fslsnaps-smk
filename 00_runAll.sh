#!/bin/bash


#default config uses hippunfoldT1    
snakemake --cores all all_centroid 
snakemake --cores all --rerun-incomplete --keep-going --notemp

for method in ashs freesurfer
do

    snakemake --cores all --configfile config/config_${method}.yml --rerun-incomplete --keep-going --notemp
done

montage results/sub-7607477/sub-7607477_hemi-R_desc-*_opacity-0_method-hippunfoldT1_snap.png results/sub-7607477/sub-7607477_hemi-[rR]*opacity-100_method-hippunfoldT1_snap.png  results/sub-7607477/sub-7607477_hemi-[rR]*opacity-100_method-ashs_snap.png results/sub-7607477/sub-7607477_hemi-[rR]*opacity-100_method-freesurfer_snap.png -geometry '1x1+0+0<' -tile 4x4 results/sub-7607477_figuremontage.png
