#!/bin/bash


#default config uses hippunfoldT1    
snakemake --cores all all_centroid 
snakemake --cores all --rerun-incomplete --keep-going --notemp

for method in ashs freesurfer
do

    snakemake --cores all --configfile config/config_${method}.yml --rerun-incomplete --keep-going --notemp
done


