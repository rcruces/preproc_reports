#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 14:13:06 2023

@author: rcruces
"""
from xhtml2pdf import pisa
import os
import argparse
import json
import glob
import nibabel as nb
import numpy as np
import matplotlib as plt
import matplotlib.pyplot as pltpy
import seaborn
from brainspace.plotting import plot_hemispheres
from brainspace.mesh.mesh_io import read_surface
from brainspace.mesh.mesh_operations import combine_surfaces
from brainspace.utils.parcellation import map_to_labels
from brainspace.vtk_interface import wrap_vtk, serial_connect
from vtk import vtkPolyDataNormals

bids='/data_/mica3/BIDS_MICs/rawdata'
sub='HC012'
ses_number='01'
sbids='sub-HC012_ses-01'
MICAPIPE='/host/yeatman/local_raid/rcruces/git_here/micapipe'

from brainspace.datasets import load_mask
mask_32k = load_mask(join=True)

c69_32k_I_lh = read_surface(MICAPIPE+'/surfaces/fsLR-32k.L.inflated.surf.gii', itype='gii')
c69_32k_I_rh = read_surface(MICAPIPE+'/surfaces/fsLR-32k.R.inflated.surf.gii', itype='gii')

def load_surface(lh, rh, with_normals=True, join=False):
    """
    Loads surfaces.

    Parameters
    ----------
    with_normals : bool, optional
        Whether to compute surface normals. Default is True.
    join : bool, optional.
        If False, return one surface for left and right hemispheres. Otherwise,
        return a single surface as a combination of both left and right.
        surfaces. Default is False.

    Returns
    -------
    surf : tuple of BSPolyData or BSPolyData.
        Surfaces for left and right hemispheres. If ``join == True``, one
        surface with both hemispheres.
    """

    surfs = [None] * 2
    for i, side in enumerate([lh, rh]):
        surfs[i] = read_surface(side)
        if with_normals:
            nf = wrap_vtk(vtkPolyDataNormals, splitting=False,
                          featureAngle=0.1)
            surfs[i] = serial_connect(surfs[i], nf)

    if join:
        return combine_surfaces(*surfs)
    return surfs[0], surfs[1]

# QC Summary
def report_qc_summary_template(jsonPath=''):
    module_description = os.path.realpath(jsonPath)
    with open( module_description ) as f:
        module_description = json.load(f)
    module = module_description["Module"]
    status = module_description["Status"]
    progress = module_description["Progress"]
    time = module_description["Processing.time"]
    threads = module_description["Threads"]
    micapipe_version = module_description["micapipeVersion"]
    date = module_description["Date"]

    # QC summary table
    report_qc_summary = (
        '<p style="font-family:Helvetica, sans-serif;font-size:12px;text-align:Left;margin-bottom:0px">'
        '<b>QC summary</b> </p>'

        '<table style="border:1px solid #666;width:100%">'
            # Status
            '<tr><td style=padding-top:4px;padding-left:3px;text-align:left><b>Status</b></td>'
            '<td style=padding-top:4px;padding-left:3px;text-align:left>{status}: {progress} steps completed</td></tr>'
            # Processing time
            '<tr><td style=padding-top:4px;padding-left:3px;text-align:left><b>Processing time</b></td>'
            '<td style=padding-top:4px;padding-left:3px;text-align:left>{time} minutes</td></tr>'
            # Thread
            '<tr><td style=padding-top:4px;padding-left:3px;text-align:left><b>Number of threads</b></td>'
            '<td style=padding-top:4px;padding-left:3px;text-align:left>{threads}</td></tr>'
            # Micapipe version
            '<tr><td style=padding-top:4px;padding-left:3px;text-align:left><b>Micapipe version</b></td>'
            '<td style=padding-top:4px;padding-left:3px;text-align:left>{micapipe_version}</td></tr>'
            # Date
            '<tr><td style=padding-top:4px;padding-left:3px;text-align:left><b>Date</b></td>'
            '<td style=padding-top:4px;padding-left:3px;text-align:left>{date}</td></tr>'
        '</table>'
    )
    return report_qc_summary.format(status=status, progress=progress, time=time, threads=threads, micapipe_version=micapipe_version, date=date)

def check_json_complete(jsonPath=''):
    module_description = os.path.realpath(jsonPath)
    with open( module_description ) as f:
        module_description = json.load(f)
    json_complete = module_description["Status"] == "COMPLETED"

    return json_complete

def report_header_template(sub='', ses_number='', dataset_name='', MICAPIPE=''):
    # Header
    report_header = (
        # Micapipe banner
        '<img id=\"top\" src=\"{MICAPIPE}/img/micapipe_long.png\" alt=\"micapipe\">'

        # Dataset name
        '<h1 style="color:#343434;font-family:Helvetica, sans-serif !important;text-align:center;margn-bottom:0">'
        '{dataset_name} <h1>'

        # Subject's ID | Session
        '<h3 style="color:#343434;font-family:Helvetica, sans-serif;text-align:center;margin-bottom:25px">'
        '<b>Subject</b>: {sub} &nbsp | &nbsp <b>Session</b>: {ses_number} </h3>'
    )
    return report_header.format(sub=sub, ses_number=ses_number, dataset_name=dataset_name, MICAPIPE=MICAPIPE)

# Header template
def report_module_header_template(module=''):
    # Module header:
    report_module_header = (
        '<p style="border:2px solid #666;padding-top:10px;padding-left:5px;background-color:#eee;font-family:Helvetica, '
        'sans-serif;font-size:14px">'
        '<b>Module: {module}</b> <p>'
    )
    return report_module_header.format(module=module)


def qc_header():
    dataset_description = os.path.realpath("%s/dataset_description.json"%(bids))
    with open( dataset_description ) as f:
        dataset_description = json.load(f)
    dataset_name = dataset_description["Name"]
    _static_block = report_header_template(sub=sub, ses_number=ses_number, dataset_name=dataset_name, MICAPIPE=MICAPIPE)

    return _static_block

subj_dir='/data/mica3/BIDS_MICs/derivatives/micapipe_v0.2.0/sub-HC012/ses-01'
swm_json=f'{subj_dir}/QC/sub-HC012_ses-01_module-SWM.json'
tmpDir="/home/bic/rcruces/Desktop/tmp_fig"

#------------------------------------------------------------------------------#
# SWM generation and mapping QC
def qc_swm(swm_json=''):
    # QC header
    _static_block = qc_header()
    _static_block +=  report_module_header_template(module='Superficial White Matter')
    
     # QC summary
    _static_block += report_qc_summary_template(swm_json)
    
    if not check_json_complete(swm_json):
        print('INCOMPLETE') # return(_static_block)
    
    # Outputs
    _static_block += (
                '<p style="font-family:Helvetica, sans-serif;font-size:12px;text-align:Left;margin-bottom:0px">'
                '<b>Main outputs</b> </p>')
    
    # SWM surfaces
    _static_block += (
                '<p style="font-family:Helvetica, sans-serif;font-size:12px;text-align:Left;margin-bottom:0px">'
                '<b>SWM surfaces</b> </p>'
                '<br />'
                '<table style="border:1px solid #666;width:100%">'
                '<tr><td style=padding-top:4px;padding-left:3px;padding-right:4px;text-align:center colspan="2"><b>fsnative</b></td></tr>')
    
    def surf_table_row(Title, png_path):
        # Add a new row to the table
        surf_row = (
        '<tr><td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:4px;text-align:center><b>"{Title}"</b></td>'
        '<td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:center><img style="display:block;width:1500px%;margin-top:0px" src="{png_path}"></td></tr>')
        return(surf_row)
    
    
    # List all the Right and Left SWM surfaces
    surf_dir = f"{subj_dir}/surf"
    swm_L_files = sorted(glob.glob(f"{surf_dir}/{sbids}_hemi-L_surf-fsnative_label-swm*gii"))
    swm_R_files = sorted(glob.glob(f"{surf_dir}/{sbids}_hemi-R_surf-fsnative_label-swm*gii"))
    
    # Load each surface and plot them and save png
    for i,_ in enumerate(swm_L_files):
        swm_L = swm_L_files[i]
        swm_R = swm_R_files[i]
        # Set the label name
        swm_label = swm_L.replace(".surf.gii","").split('label-')[1]
        # Load the SWM surface
        lhSWM, rhSWM = load_surface(swm_L, swm_R, with_normals=True, join=False)
        
        # Plot the surfaces SWM
        #dsize = (900, 250)
        #display = Display(visible=0, size=dsize)
        #display.start()
        swm_png=f"{tmpDir}/{sbids}_space-nativepro_surf-fsnative_label-{swm_label}.png"
        plot_hemispheres(lhSWM, rhSWM, size=(900, 250), zoom=1.25, embed_nb=True, interactive=False, share='both',
                         nan_color=(0, 0, 0, 1), color_range=(-1,1), transparent_bg=False,
                         screenshot = True, offscreen=True, filename = swm_png)
        _static_block += surf_table_row(swm_label, swm_png)
        #display.stop()
    
    _static_block += ('</table>')
        
    #------------------------------------------------------------------------------#
    # SWM maps on surfs
    _static_block += (
                '<p style="font-family:Helvetica, sans-serif;font-size:12px;text-align:Left;margin-bottom:0px">'
                '<b>SWM maps</b> </p>')
    
    # List all the Right and Left SWM surfaces
    map_dir = f"{subj_dir}/maps"
    swm_files = sorted(glob.glob(f"{map_dir}/{sbids}_hemi-L_surf-fsLR-32k_label-swm*.func.gii"))
    
    # Get the unique maps IDs
    maps_str = list(set([file.split('mm_')[1][:-9] for file in swm_files]))
    
    # Get the unique surfaces
    surf_str = sorted(list(set([file.split('label-')[1].split('_')[0] for file in swm_files])))
    
    for measure in maps_str:
            swm_surf_table = '<br>'
            swm_surf_table += (
                '<table style="border:1px solid #666;width:100%">'
                     '<tr><td style=padding-top:4px;padding-left:3px;padding-right:4px;text-align:center colspan="2"><b>{measure} (fsLR-32k)</b></td></tr>'
             ).format(measure=measure)
    
            for i, surface in enumerate(surf_str):
                measure_c69_32k_lh = f"{map_dir}/{sbids}_hemi-L_surf-fsLR-32k_label-{surface}_{measure}.func.gii"
                measure_c69_32k_rh = f"{map_dir}/{sbids}_hemi-R_surf-fsLR-32k_label-{surface}_{measure}.func.gii"
                measure_c69_32k_png = f"{tmpDir}/{sbids}_surf-fsLR-32k_label-{surface}_{measure}.png"
                f = np.concatenate((nb.load(measure_c69_32k_lh).darrays[0].data, nb.load(measure_c69_32k_rh).darrays[0].data), axis=0)
                
                if i == 0: measure_crange=(np.quantile(f[mask_32k], 0.01), np.quantile(f[mask_32k], 0.98))
                #display = Display(visible=0, size=(900, 250))
                #display.start()
                # Replace values in f with NaN where mask_32k is False
                f[mask_32k == False] = np.nan
                # PLot the values
                plot_hemispheres(c69_32k_I_lh, c69_32k_I_rh, array_name=f, size=(900, 250), color_bar='bottom', zoom=1.25, embed_nb=True, interactive=False, share='both',
                                 nan_color=(0, 0, 0, 1), color_range=measure_crange, cmap='mako', transparent_bg=False,
                                 screenshot = True, offscreen=True, filename = measure_c69_32k_png)
                #display.stop()
                swm_surf_table += (
                         '<tr><td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:4px;text-align:center><b>{surface}</b></td>'
                         '<td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:center><img style="display:block;width:1500px%;margin-top:0px" src="{measure_c69_32k_png}"></td></tr>'
                ).format(surface=surface,
                    measure_c69_32k_png=measure_c69_32k_png
                )
            _static_block += swm_surf_table
    
            _static_block += '</table>'
            
    return _static_block

qc_swm(swm_json)
