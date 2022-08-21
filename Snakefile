import numpy as np
from snakebids import bids
from snakebids.utils.snakemake_io import glob_wildcards
from snakemake.io import get_wildcard_names


configfile: 'config_hippunfold.yml'

#a kind of lightweight snakebids below:


wildcard_names = list(get_wildcard_names(config['seg_path']))

#with glob_wildcards it seems we may need to always hard-code the choice of wildcards.. unless either snakebids is used, or modify glob_wildcards to return a dict..
segnamedtuple = glob_wildcards(config['seg_path'])

seg_zip_list = dict()
seg_wildcards = dict()

#populate the seg_zip_list and seg_wildcards:
for name in wildcard_names:
    seg_zip_list[name] = getattr(segnamedtuple, name)
    seg_wildcards[name] = f'{{{name}}}'



def get_coords(wildcards, input):

    centroid = np.loadtxt(input.centroid)    
    x = float(wildcards.xoffset) + centroid[0]
    y = float(wildcards.yoffset) + centroid[1]
    z = float(wildcards.zoffset) + centroid[2]
    return f'{x} {y} {z}'



rule all:
    input:
        expand(
        expand(
            bids(root='results',
                suffix='opacitymontage.png',
                desc='{desc}',
                **seg_wildcards),
            zip,
            **seg_zip_list,
            allow_missing=True),
                desc=config['slice_montages'].keys(),
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
            suffix='slicemontage.png',
            desc='{desc}',
            opacity='{opacity}',
            **seg_wildcards),
                opacity=['0','100'],allow_missing=True)
    params:
        tile='1x2',
        geometry="'1x1+0+0<'"
    output:
            bids(root='results',
            suffix='opacitymontage.png',
            desc='{desc}',
            **seg_wildcards),
    shell:
        'montage {input} -geometry {params.geometry} -tile {params.tile} {output}'

        
def get_input_slices(wildcards):

    offset = dict()
    for ax in ['x','y','z']:
    
        start = config['slice_montages'][wildcards.desc][ax]['start']
        stop = config['slice_montages'][wildcards.desc][ax]['stop']
        slices = config['slice_montages'][wildcards.desc][ax]['slices']
        offset[ax] = [f'{num}' for num in np.linspace(start,stop,slices)]

    return expand(
                bids(root='results',suffix='snap.png',
                x='{x}',
                y='{y}',
                z='{z}',
                desc='{desc}',
                opacity='{opacity}',
                **seg_wildcards),
                    x=offset['x'],
                    y=offset['y'],
                    z=offset['z'],
                    allow_missing=True)


 
rule montage_slices:
    input:
        slices = get_input_slices
    params:
        tile=lambda wildcards, input: '{N}x1'.format(N=len(input)),
        geometry='800x600'
    output:
        bids(root='results',
            suffix='slicemontage.png',
            desc='{desc}',
            opacity='{opacity}',
            **seg_wildcards)
    shell:
        'montage {input} -geometry {params.geometry} -tile {params.tile} {output}'

rule get_slice_centroid:
    input:
        seg = config['seg_path'].format(**seg_wildcards),
    output:
        centroid = bids(root='results',suffix='centroid.txt',**seg_wildcards),
    shell:
        'fslstats  {input} -c > {output}'


def get_hide_slices(wildcards):
    if wildcards.desc == 'coronal':
        return '--hidex --hidez'
    if wildcards.desc == 'axial':
        return '--hidex --hidey'
    if wildcards.desc == 'sagittal':
        return '--hidey --hidez'



        
    

rule gen_snap:
    """ generates a snapshot with location relative to centroid of the segmentation """
    input:
        mri = config['mri_path'].format(**seg_wildcards),
        seg = config['seg_path'].format(**seg_wildcards),
        centroid = bids(root='results',suffix='centroid.txt',**seg_wildcards),
    params:
        coords = get_coords,
        label_opacity = '{opacity}',
        lut = config['lut'],
        hide_slices = get_hide_slices,
        zoom = lambda wildcards: config['slice_montages'][wildcards.desc]['zoom']
    output:
        bids(root='results',suffix='snap.png',
                x='{xoffset}',
                y='{yoffset}',
                z='{zoffset}',
                desc='{desc}',
                opacity='{opacity}',
                **seg_wildcards)
    shell:
        "fsleyes render -of {output}"
        " --scene ortho"
        " --worldLoc {params.coords}"
        " --displaySpace {input.mri}"
        " --xzoom {params.zoom}"
        " --yzoom {params.zoom}"
        " --zzoom {params.zoom}"
        " --layout horizontal"
        " {params.hide_slices}"
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
        " --lut {params.lut}" #random"
        " --outlineWidth 0"
        " --volume 0"



