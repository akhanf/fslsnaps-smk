import numpy as np
from snakebids import bids
from snakebids.utils.snakemake_io import glob_wildcards
from snakemake.io import get_wildcard_names



seg_path = '../hippunfold_highresT2/hippunfold/sub-{subject}/anat/sub-{subject}_hemi-{hemi}_space-cropT2w_desc-subfields_atlas-{atlas}_dseg.nii.gz' 

#note: the wildcards in mri_path must be a subset of those in seg_path 
mri_path = '../hippunfold_highresT2/hippunfold/sub-{subject}/anat/sub-{subject}_desc-preproc_T2w.nii.gz'



#a kind of lightweight snakebids below:

wildcard_names = list(get_wildcard_names(seg_path))

#with glob_wildcards it seems we may need to always hard-code the choice of wildcards.. unless either snakebids is used, or modify glob_wildcards to return a dict..
segnamedtuple = glob_wildcards(seg_path)

seg_zip_list = dict()
seg_wildcards = dict()

#populate the seg_zip_list and seg_wildcards:
for name in wildcard_names:
    seg_zip_list[name] = getattr(segnamedtuple, name)
    seg_wildcards[name] = f'{{{name}}}'


print(seg_zip_list)
print(seg_wildcards)

def get_coords(wildcards, input):

    centroid = np.loadtxt(input.centroid)    
    x = float(wildcards.xoffset) + centroid[0]
    y = float(wildcards.yoffset) + centroid[1]
    z = float(wildcards.zoffset) + centroid[2]
    return f'{x} {y} {z}'


def get_input_slices(wildcards):

    start = float(wildcards.start)
    stop = float(wildcards.stop)
    slices = int(wildcards.slices)
    yoffset=[f'{num}' for num in np.linspace(start,stop,slices)]

    return expand(
                bids(root='results',suffix='snap.png',
                x='{xoffset}',
                y='{yoffset}',
                z='{zoffset}',
                opacity='{opacity}',
                **seg_wildcards),
                    xoffset='0',yoffset=yoffset,zoffset='0',
                    allow_missing=True)


rule all:
    input:
        expand(
        expand(
            bids(root='results',
                suffix='montage.png',
                desc='{desc}',
                start='{start}',
                stop='{stop}',
                slices='{slices}',
                **seg_wildcards),
            zip,
            **seg_zip_list,
            allow_missing=True),
                desc='onoff',
                start='-15',
                stop='15',
                slices=5
)

rule all_centroids:
    input:
        expand(bids(root='results',suffix='centroid.txt',**seg_wildcards),
            zip,
            **seg_zip_list,
            allow_missing=True),


#TODO: could easily make an animated rule (from opacity 0 to 100 to 0)

rule montage_seg_with_mri:
    input:
        expand(
            bids(root='results',
            suffix='montage.png',
            opacity='{opacity}',
            start='{start}',
            stop='{stop}',
            slices='{slices}',
            **seg_wildcards),
                opacity=['0','100'],allow_missing=True)
    params:
        tile='1x2',
        geometry=lambda wildcards: '{xdim}x600'.format(xdim=800*int(wildcards.slices))
    output:
            bids(root='results',
            suffix='montage.png',
            desc='onoff',
            start='{start}',
            stop='{stop}',
            slices='{slices}',
            **seg_wildcards),
    shell:
        'montage {input} -geometry {params.geometry} -tile {params.tile} {output}'

        

 
rule montage_coronals:
    input:
        slices = get_input_slices 
    params:
        tile=lambda wildcards, input: '{N}x1'.format(N=len(input)),
        geometry='800x600'
    output:
        bids(root='results',
            suffix='montage.png',
            opacity='{opacity}',
            start='{start}',
            stop='{stop}',
            slices='{slices}',
            **seg_wildcards)
    shell:
        'montage {input} -geometry {params.geometry} -tile {params.tile} {output}'

rule get_slice_centroid:
    input:
        seg = seg_path.format(**seg_wildcards),
    output:
        centroid = bids(root='results',suffix='centroid.txt',**seg_wildcards),
    shell:
        'fslstats  {input} -c > {output}'


rule gen_snap:
    """ generates a snapshot with location relative to centroid of the segmentation """
    input:
        mri = mri_path.format(**seg_wildcards),
        seg = seg_path.format(**seg_wildcards),
        centroid = bids(root='results',suffix='centroid.txt',**seg_wildcards),
    params:
        coords = get_coords,
        label_opacity = '{opacity}'
    output:
        bids(root='results',suffix='snap.png',
                x='{xoffset}',
                y='{yoffset}',
                z='{zoffset}',
                opacity='{opacity}',
                **seg_wildcards)
    shell:
        "fsleyes render -of {output}"
        " --scene ortho"
        " --worldLoc {params.coords}"
        " --displaySpace {input.mri}"
        " --xzoom 2386.6666666666665"
        " --yzoom 2386.6666666666665"
        " --zzoom 2386.6666666666665"
        " --layout horizontal"
        " --hidex"
        " --hidez"
        " --hideCursor"
        " --hideLabels"
        " --bgColour 0.0 0.0 0.0"
        " --fgColour 1.0 1.0 1.0"
        " --cursorColour 0.0 1.0 0.0"
        " --colourBarLocation top"
        " --colourBarLabelSide top-left"
        " --colourBarSize 100.0"
        " --labelSize 12"
        " --performance 3"
        " --movieSync {input.mri}"
        " --name 'mri'"
        " --overlayType volume"
        " --alpha 100.0"
        " --brightness 61.692289455328854"
        " --contrast 73.38457891065771"
        " --cmap greyscale"
        " --negativeCmap greyscale"
        " --displayRange 0.0 700.0"
        " --clippingRange 0.0 1328.1773706054687"
        " --modulateRange 0.0 1315.027099609375"
        " --gamma 0.0"
        " --cmapResolution 256"
        " --interpolation spline"
        " --interpolateCmaps"
        " --numSteps 100"
        " --blendFactor 0.1"
        " --smoothing 0"
        " --resolution 100"
        " --numInnerSteps 10"
        " --clipMode intersection"
        " --volume 0 {input.seg}"
        " --name 'seg'"
        " --overlayType label"
        " --alpha {params.label_opacity}"
        " --brightness 49.75000000000001"
        " --contrast 49.90029860765409"
        " --lut random"
        " --outlineWidth 0"
        " --volume 0"



