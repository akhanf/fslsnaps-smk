#!/bin/bash

for method in hippunfoldt1 ashs freesurfer
do

    snakemake --cores all --configfile config/config_${method}.yml 
done
