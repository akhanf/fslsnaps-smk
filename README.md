# fslsnaps-smk


1. Edit these specific fields in the config.yml:

  - `method:` name of the segmentation method to include in the final file names
  - `lut:` path to look-up table to use for coloring the segmentation
  - `custom_path:` path (with wildcards) to use for the mri and seg. This appears in the `mri` and the `seg` section.


2. Run the `all_centroids` target rule first to collect the centroid and slice information for all the segmentations.
    
```
snakemake all_centroids --cores all
```

3. Run the workflow with the default target rule (`all`) to generate all the snapshots for the segmentations (final output is the flipbook pdf)

```
snakemake all --cores all
```
which is the same as:
```
snakemake --cores all
```

4. If you want to generate flipbooks for additional segmentations, create another config file, and customize it as in Step 1, then run 
snakemake with the `--configfile PATH_TO_CONFIG_FILE` option. Be sure to change the `method` and seg `custom_path` field. Note that we skip Step 2 to ensure identical slices are used. E.g. if the additional config file is placed in `config/config_freesurfer.yml`:

```
snakemake --configfile config/config_freesurfer.yml --cores all
```


