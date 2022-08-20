import numpy as np

# to do:
#  - use config file for other optoins (also maybe start stop, nslices too)
#  - figure out best way to iterate over subject images


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

    return expand('test_x0_y{y}_z0_opacity-{{opacity}}.png',y=yoffset)

 
rule montage_coronals:
    input:
        slices = get_input_slices 
    params:
        tile=lambda wildcards, input: '{N}x1'.format(N=len(input)),
        geometry='800x600'
    output:
        'test_montage_coronal_linspace_{start}_{stop}_{slices}_opacity-{opacity}.png'
    shell:
        'montage {input} -geometry {params.geometry} -tile {params.tile} {output}'

rule get_slice_centroid:
    input:
        '/local/scratch/hippunfold-fs-ashs-comparison/hippunfold_highresT2/hippunfold/sub-9992517/anat/sub-9992517_hemi-L_space-cropT2w_desc-subfields_atlas-bigbrain_dseg.nii.gz'
    output:
        'centroid.txt'
    shell:
        'fslstats  /local/scratch/hippunfold-fs-ashs-comparison/hippunfold_highresT2/hippunfold/sub-9992517/anat/sub-9992517_hemi-L_space-cropT2w_desc-subfields_atlas-bigbrain_dseg.nii.gz -c > {output}'


rule gen_snap:
    """ generates a snapshot with location relative to centroid of the segmentation """
    input:
        mri = '/local/scratch/hippunfold-fs-ashs-comparison/hippunfold_highresT2/hippunfold/sub-9992517/anat/sub-9992517_hemi-L_space-cropT2w_desc-subfields_atlas-bigbrain_dseg.nii.gz',
        " --displaySpace /local/scratch/hippunfold-fs-ashs-comparison/hippunfold_highresT2/hippunfold/sub-9992517/anat/sub-9992517_desc-preproc_T2w.nii.gz"
        seg = '/local/scratch/hippunfold-fs-ashs-comparison/hippunfold_highresT2/hippunfold/sub-9992517/anat/sub-9992517_hemi-L_space-cropT2w_desc-subfields_atlas-bigbrain_dseg.nii.gz',
        centroid = 'centroid.txt'
    params:
        coords = get_coords,
        label_opacity = '{opacity}'
    output:
        'test_x{xoffset}_y{yoffset}_z{zoffset}_opacity-{opacity}.png'
    shell:
        "fsleyes render -of {output}"
        " --scene ortho"
        " --worldLoc {params.coords}"
        " --displaySpace {input.mri}"
        " --xcentre  0.03418 -0.22980"
        " --ycentre  0.35016 -0.26185"
        " --zcentre  0.31203  0.03418"
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



