import numpy as np
from snakebids import bids
from snakebids.utils.snakemake_io import glob_wildcards
from snakemake.io import get_wildcard_names
from snakebids import filter_list
from snakebids import generate_inputs
from snakebids import get_wildcard_constraints

configfile: 'config/config.yml'

inputs = generate_inputs(bids_dir=config['bids_dir'],
                                    pybids_inputs=config['pybids_inputs'],
                                    use_bids_inputs=True)



#add in wildcard constraints
wildcard_constraints: **get_wildcard_constraints(config['pybids_inputs'])

wildcard_constraints:
    method='[a-zA-Z0-9]+',


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
        expand(bids(root='results',
                suffix='flipbook.pdf',
                hemi='{hemi}',
                method='{method}',
                include_subject_dir=False),
            hemi=inputs['seg'].input_lists['hemi'],
            method=config['method'])

rule all_centroid:
    input:
        expand(bids(root='results',suffix='centroid.txt',**inputs['seg'].input_wildcards),
            zip,
            **inputs['seg'].input_zip_lists,
            allow_missing=True),
        expand(bids(root='results',suffix='vox2ras.txt',**inputs['seg'].input_wildcards),
            zip,
            **inputs['seg'].input_zip_lists,
            allow_missing=True),


#combine into pdfs for flipping through -- same hemi and method 

def get_inputs_create_pdf(wildcards):
    
    pngs = expand(
                    bids(root='results',
                        suffix='viewmontage.png',
                        method='{method}'.format(method=wildcards.method),
                        **inputs['seg'].input_wildcards),
                    zip,
                    **filter_list(inputs['seg'].input_zip_lists,wildcards)
                    )
    return pngs
            


rule create_pdf:
    input:
        get_inputs_create_pdf
    output:
        bids(root='results',
                suffix='flipbook.pdf',
                hemi='{hemi}',
                method='{method}',
                include_subject_dir=False)
    shell:
        'convert {input} {output}'




#TODO: could easily make an animated rule (from opacity 0 to 100 to 0)





rule montage_views:
    input:
        expand(bids(root='results',
            suffix='opacitymontage.png',
            desc='{desc}',
            method='{method}',
            **inputs['seg'].input_wildcards),
                desc=config['slice_montages'].keys(),
                allow_missing=True)
    params:
        tile=lambda wildcards, input: '1x{N}'.format(N=len(input)),
        geometry="'1x1+0+0<'" 
    output:
        bids(root='results',
            suffix='viewmontage.png',
            method='{method}',
            **inputs['seg'].input_wildcards),
    shadow: 'minimal'
    shell:
        "montage {input}  -geometry {params.geometry} -tile {params.tile} temp.png && " 
        "convert temp.png -background black -gravity North -splice 0x40 -fill white -pointsize 30 -annotate +0+2 '{wildcards}'"
         " {output}"
   
    

rule montage_seg_with_mri:
    """ stack up slice with no seg (opacity 0) over with seg (opacity 100)"""
    input:
        expand(
            bids(root='results',
            suffix='slicemontage.png',
            desc='{desc}',
            opacity='{opacity}',
            method='{method}',
            **inputs['seg'].input_wildcards),
                opacity=['0',config['opacity']],allow_missing=True)
    params:
        tile=lambda wildcards, input: '1x{N}'.format(N=len(input)),
        geometry="'1x1+0+0<'"
    output:
            bids(root='results',
            suffix='opacitymontage.png',
            desc='{desc}',
            method='{method}',
            **inputs['seg'].input_wildcards),
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
                method='{method}',
                **inputs['seg'].input_wildcards),
                    x=offset['x'],
                    y=offset['y'],
                    z=offset['z'],
                    allow_missing=True)


 
rule montage_slices:
    input:
        slices = get_input_slices
    params:
        tile=lambda wildcards, input: '{N}x1'.format(N=len(input)),
        geometry="'1x1+0+0<'" 
    output:
        bids(root='results',
            suffix='slicemontage.png',
            desc='{desc}',
            opacity='{opacity}',
            method='{method}',
            **inputs['seg'].input_wildcards)
    shell:
        'montage {input} -geometry {params.geometry} -tile {params.tile} {output}'


hemi_standardize = {'left': 'L', 'right': 'R', 'lh':'L', 'rh':'R','L': 'L', 'R':'R'}

def get_seg_ref(wildcards):
    hemi = hemi_standardize[wildcards.hemi]
    return config['seg_ref_path'].format(hemi=hemi,subject=wildcards.subject)

rule get_slice_centroid:
    input:
        seg = get_seg_ref,
    output:
        centroid = bids(root='results',suffix='centroid.txt',**inputs['seg'].input_wildcards)
    shell:
        'fslstats  {input} -c > {output}'

rule get_vox_to_ras:
    input:
        seg = get_seg_ref,
    output:
        vox2ras = bids(root='results',suffix='vox2ras.txt',**inputs['seg'].input_wildcards),
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
    return bids(root='results',suffix='centroid.txt',**inputs['seg'].input_wildcards).format(hemi=hemi,subject=wildcards.subject)


def get_vox2ras_txt(wildcards):
    hemi = hemi_standardize[wildcards.hemi]
    return bids(root='results',suffix='vox2ras.txt',**inputs['seg'].input_wildcards).format(hemi=hemi,subject=wildcards.subject)

       
    

rule gen_snap:
    """ generates a snapshot with location relative to centroid of the segmentation """
    input:
        mri = inputs['mri'].input_path.format(**inputs['seg'].input_wildcards),
        seg = inputs['seg'].input_path.format(**inputs['seg'].input_wildcards),
        centroid = get_centroid_txt,
        vox2ras = get_vox2ras_txt,
        lut = config['lut'],
    params:
        img_size = config['size'],
        coords = get_coords,
        label_opacity = '{opacity}',
        displayrange = config['displayrange_percentiles'],
        hide_slices = get_hide_slices,
        zoom = lambda wildcards: config['slice_montages'][wildcards.desc]['zoom']
    output:
        temp(bids(root='results',suffix='snap.png',
                x='{xoffset}',
                y='{yoffset}',
                z='{zoffset}',
                desc='{desc}',
                opacity='{opacity}',
                method='{method}',
                **inputs['seg'].input_wildcards))
    shell:
        "xvfb-run -a fsleyes render -of {output}"
        " --size {params.img_size}"
        " --initialDisplayRange {params.displayrange}"
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
        " --cmap greyscale"
        " --negativeCmap greyscale"
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
        " --lut {input.lut}" 
        " --outlineWidth 0"
        " --volume 0"



