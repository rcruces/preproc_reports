#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 14:36:26 2023

@author: rcruces
"""

import os
import glob
import nibabel as nib

subj_func='/data_/mica3/BIDS_PNI/derivatives/micapipe_v0.2.0/'

gii = glob.glob(subj_func+'sub-*/ses-*/func/desc-me_task-*/volumetric/*SNR*.shape.gii')

for gii_file in gii:
    # Split the string by '/'
    parts = gii_file.split('/')[11].replace('_space-func_desc-me_tSNR.shape.gii','')
    task = gii_file.split('/')[9].replace('desc-me_task-','')
    data = nib.load(gii_file).darrays[0].data
    print(f"{data.shape} : {parts} {task}")

from brainspace.plotting import plot_hemispheres
from brainspace.mesh.mesh_io import read_surface
from brainspace.vtk_interface import wrap_vtk, serial_connect
from vtk import vtkPolyDataNormals
from brainspace.mesh.mesh_operations import combine_surfaces
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


# Nativepro fsnative white
sub="sub-HC062"
ses="ses-02"
surf_id=f"{sub}_{ses}"
surf=f"/data/mica3/BIDS_MICs/derivatives/micapipe_v0.2.0/{sub}/{ses}/surf"

# Surfaces
def surf_points(label, Surf=False):
    if Surf != False:
        Nom = f"surf-fsnative_label-{label}.surf.gii"
        orig="surf-fsnative"
    else:
        Nom = f"space-nativepro_surf-fsnative_label-{label}.surf.gii"
        orig="space-nativepro"
        
    lh, rh = load_surface(f"{surf}/{surf_id}_hemi-L_{Nom}" ,
                                f"{surf}/{surf_id}_hemi-R_{Nom}" , with_normals=True, join=False)
    lh_points = str(lh.points.shape[0])
    rh_points = str(rh.points.shape[0])
    print(f"{lh_points}, {rh_points} : {orig}-{label}")
        
# Number of points per subjects
surf_points("white")
surf_points("pial")
surf_points("midthickness")
surf_points("midthickness", "surf")
surf_points("sphere", "surf")

# change working diretory to 
derivatives="/data/mica3/BIDS_MICs/derivatives/micapipe_v0.2.0"
os.chdir(derivatives)

# List all the directories (ls -d sub*/ses*)
import glob
files=sorted(glob.glob("sub-HC062/ses-*/QC/*_module-post_structural.json"))

# For loop over each subject:
for file in files:
    sub=file.split('/')[0]
    ses=file.split('/')[1]
    surf_id=f"{sub}_{ses}"
    surf=f"/data/mica3/BIDS_MICs/derivatives/micapipe_v0.2.0/{sub}/{ses}/surf"
    print("-------------------------------------------------------------------")
    print(f"{surf_id}")
    # Number of points per subjects
    surf_points("white")
    surf_points("pial")
    surf_points("midthickness")
    surf_points("midthickness", "surf")
    surf_points("sphere", "surf")
    
# Plot surfaces
lhWc, rhWc = load_surface(f"{surf}/{surf_id}_hemi-L_surf-fsnative_label-swm03.surf.gii" ,
                            f"{surf}/{surf_id}_hemi-R_surf-fsnative_label-swm03.surf.gii" , with_normals=True, join=False)

plot_hemispheres(lhWc, rhWc,  size=(900, 250), zoom=1.25, embed_nb=True, interactive=False, share='both',
                     nan_color=(0, 0, 0, 1), color_range=(-1,1), cmap='Greys', transparent_bg=False,
                     screenshot = False)