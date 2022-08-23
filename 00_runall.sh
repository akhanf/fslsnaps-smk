#!/bin/bash
snakemake --cores all all_centroid && for method in ashs hippunfold freesurfer hippunfoldT1; do snakemake --config method=$method --cores all --keep-going --rerun-incomplete; done
