# fslsnaps-smk


1. Edit the `seg` dict in the config.yml to point to your segmentations specific fields in the config.yml:

  - `method:` keys in `seg` are the names of the segmentation method, to include in the final file names
  - `lut:` path to look-up table to use for coloring the segmentation
  - `path:` path (with wildcards) to use for the mri and seg. This appears in the `mri` and the `seg` section.


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


