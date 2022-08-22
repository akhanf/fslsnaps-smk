import numpy as np
from snakebids import bids
from snakebids.utils.snakemake_io import glob_wildcards
from snakemake.io import get_wildcard_names


configfile: 'config_hippunfold.yml'

#a kind of lightweight snakebids below:
wildcard_constraints:
    method='[a-zA-Z0-9]+',
    hemi='[a-zA-Z0-9]+',
    subject='[a-zA-Z0-9]+'

def get_zip_list_and_wildcards(seg_path):
    wildcard_names = list(get_wildcard_names(seg_path))

    #with glob_wildcards it seems we may need to always hard-code the choice of wildcards.. unless either snakebids is used, or modify glob_wildcards to return a dict..
    segnamedtuple = glob_wildcards(seg_path)

    seg_zip_list = dict()
    seg_wildcards = dict()

    #populate the seg_zip_list and seg_wildcards:
    for name in wildcard_names:
        seg_zip_list[name] = getattr(segnamedtuple, name)
        seg_wildcards[name] = f'{{{name}}}'

    return (seg_zip_list, seg_wildcards)


(seg_zip_list,seg_wildcards) = get_zip_list_and_wildcards(config['seg_path'])
(seg_ref_zip_list,seg_ref_wildcards) = get_zip_list_and_wildcards(config['seg_ref_path'])



def get_coords(wildcards, input):

    centroid = np.loadtxt(input.centroid).reshape([3,1])    
    vox2ras = np.loadtxt(input.vox2ras).reshape([4,4])

    #offset from wildcards is in vox space
    #centroid we have is in ras space
    #so we would need to:
    #  1. bring centroid to vox space
    #  2. add the offset
    #  3. bring it back to ras space

    #  if we had voxel-space coords, then we would:
    #  1. add the offset
    #  2. bring it back to ras space

    # let's stick with centroids in ras space for now 

    centroid_vec = np.vstack((centroid,np.array([1.0]))) # concat 1 to make homog 4x1 vec
#    print(centroid_vec)
    centroid_vec_voxspace = np.linalg.inv(vox2ras) @ centroid_vec
#    print(centroid_vec_voxspace)
    offset_vec = np.vstack((float(wildcards.xoffset),float(wildcards.yoffset),float(wildcards.zoffset),1.0))
#    print(offset_vec)
    coords_voxspace = centroid_vec_voxspace.copy()
    coords_voxspace[:3,0] = centroid_vec_voxspace[:3,0] + offset_vec[:3,0] #skip the homog coord
#    print(coords_voxspace)
    coords_ras = np.squeeze(vox2ras @ coords_voxspace)
#    print(coords_ras)
    coords_string = f'{coords_ras[0]} {coords_ras[1]} {coords_ras[2]}' 
    return coords_string



rule all:
    input:
        expand(
            bids(root='results',
                suffix='viewmontage.png',
                **seg_wildcards),
            zip,
            **seg_zip_list
        )

rule all_centroid:
    input:
        expand(bids(root='results',suffix='centroid.txt',**seg_ref_wildcards),
            zip,
            **seg_ref_zip_list,
            allow_missing=True),
        expand(bids(root='results',suffix='vox2ras.txt',**seg_ref_wildcards),
            zip,
            **seg_ref_zip_list,
            allow_missing=True),


#TODO: could easily make an animated rule (from opacity 0 to 100 to 0)



#output rule needs seg_wildcards with hemi removed
#seg_wildcards_nohemi = seg_wildcards.copy()
#del seg_wildcards_nohemi['hemi']

#rule montage_hemis:
#    """ stack these left to right. removes the hemi wildcard -- doesn't look great though.."""
#
#    input: 
#        expand(bids(root='results',
#            suffix='viewmontage.png',
#            **seg_wildcards),
#            hemi=['R','L'],allow_missing=True)
#    output:
#        bids(root='results',
#            suffix='hemimontage.png',
#            **seg_wildcards_nohemi)
#    params:
#        tile='2x1',
#        geometry="'1x1+0+0<'" 
#    shell:
#        'montage {input} -geometry {params.geometry} -tile {params.tile} {output}'
   
 




rule montage_views:
    input:
        expand(bids(root='results',
            suffix='opacitymontage.png',
            desc='{desc}',
            **seg_wildcards),
                desc=config['slice_montages'].keys(),
                allow_missing=True)
    params:
        tile='1x2',
        geometry="'1x1+0+0<'" 
    output:
        bids(root='results',
            suffix='viewmontage.png',
            **seg_wildcards),
    shell:
        'montage {input} -geometry {params.geometry} -tile {params.tile} {output}'
   
    

rule montage_seg_with_mri:
    """ stack up slice with no seg (opacity 0) over with seg (opacity 100)"""
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


hemi_standardize = {'left': 'L', 'right': 'R', 'L': 'L', 'R':'R'}

def get_seg_ref(wildcards):
    hemi = hemi_standardize[wildcards.hemi]
    return config['seg_ref_path'].format(hemi=hemi,subject=wildcards.subject)

rule get_slice_centroid:
    input:
        seg = get_seg_ref,
    output:
        centroid = bids(root='results',suffix='centroid.txt',**seg_ref_wildcards)
    shell:
        'fslstats  {input} -c > {output}'

rule get_vox_to_ras:
    input:
        seg = get_seg_ref,
    output:
        vox2ras = bids(root='results',suffix='vox2ras.txt',**seg_ref_wildcards),
    shell:
        'fslorient -getqform {input} > {output}'



def get_hide_slices(wildcards):
    if wildcards.desc == 'coronal':
        return '--hidex --hidez'
    if wildcards.desc == 'axial':
        return '--hidex --hidey'
    if wildcards.desc == 'sagittal':
        return '--hidey --hidez'


def get_centroid_txt(wildcards):
    hemi = hemi_standardize[wildcards.hemi]
    return bids(root='results',suffix='centroid.txt',**seg_ref_wildcards).format(hemi=hemi,subject=wildcards.subject)


def get_vox2ras_txt(wildcards):
    hemi = hemi_standardize[wildcards.hemi]
    return bids(root='results',suffix='vox2ras.txt',**seg_ref_wildcards).format(hemi=hemi,subject=wildcards.subject)

       
    

rule gen_snap:
    """ generates a snapshot with location relative to centroid of the segmentation """
    input:
        mri = config['mri_path'].format(**seg_wildcards),
        seg = config['seg_path'].format(**seg_wildcards),
        centroid = get_centroid_txt,
        vox2ras = get_vox2ras_txt
    params:
        coords = get_coords,
        label_opacity = '{opacity}',
        lut = config['lut'],
        hide_slices = get_hide_slices,
        zoom = lambda wildcards: config['slice_montages'][wildcards.desc]['zoom']
    output:
        temp(bids(root='results',suffix='snap.png',
                x='{xoffset}',
                y='{yoffset}',
                z='{zoffset}',
                desc='{desc}',
                opacity='{opacity}',
                **seg_wildcards))
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



