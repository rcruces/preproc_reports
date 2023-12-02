#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 14:39:01 2023

@author: rcruces
"""
import glob
import nibabel as nb
import numpy as np
from brainspace.mesh.mesh_io import read_surface
from brainspace.gradient import GradientMaps
from brainspace.plotting import plot_hemispheres

# micapipe directory
micapipe='/data_/mica1/01_programs/micapipe-v0.2.0'
# Shape of the fsLR-5k matrices
N5k = 9684

# -----------------------------------------------------------------------------
# Mask medial wall mask
# fsLR-5k mask
mask_lh = nb.load(micapipe + '/surfaces/fsLR-5k.L.mask.shape.gii').darrays[0].data
mask_rh = nb.load(micapipe + '/surfaces/fsLR-5k.R.mask.shape.gii').darrays[0].data
mask_5k = np.concatenate((mask_lh, mask_rh), axis=0)

# -----------------------------------------------------------------------------
# Surfaces to plot
# Load fsLR-5k inflated surface
i5_lh = read_surface(micapipe + '/surfaces/fsLR-5k.L.inflated.surf.gii', itype='gii')
i5_rh = read_surface(micapipe + '/surfaces/fsLR-5k.R.inflated.surf.gii', itype='gii')

# -----------------------------------------------------------------------------
# FILES of SC from PNI 7T all subjects
# List all the files
# sc_file = sorted(glob.glob(f"sub-*/ses-*/dwi/connectomes/*_surf-fsLR-5k_desc-iFOD2-40M-SIFT2_full-connectome.shape.gii"))
# #ed_file = sorted(glob.glob(f"sub-*/ses-*/dwi/connectomes/*_surf-fsLR-5k_desc-iFOD2-40M-SIFT2_full-edgeLengths.shape.gii"))

# Load the files
sc_file = np.load('/host/yeatman/local_raid/rcruces/data/tmp/sc_file.npy')

# -----------------------------------------------------------------------------
# All SC from PNI 7T all subjects 
# Data was loaded like this:
def load_sc(File):
    """Loads and process a structura connectome"""
    
    # load the matrix
    mtx_sc = nb.load(File).darrays[0].data
    
    # Mirror the matrix & log
    mtx_sc = np.triu(mtx_sc,1)+mtx_sc.T
    
    # Replace NaN with 0
    mtx_sc[np.isnan(mtx_sc)] = 0
    
    # Replace negative infinite with 0
    mtx_sc[np.isneginf(mtx_sc)] = 0
    
    return mtx_sc

# Loads all the MPC fsLR-5k matrices
sc_5k=np.empty([N5k, N5k, len(sc_file)], dtype=float)

# Load the SC and weight it by the edge length
for i, f in enumerate(sc_file):
    print(f"{i} - {f}")
    #sc_5k[:,:,i] = load_sc(f) / load_sc(ed_file[i]) # Weighted SC by the Edgelengths
    sc_5k[:,:,i] = load_sc(f)

# Load ALL the connectomes (25G!!!! 33 matrices)
# SC_5k = np.load('/host/yeatman/local_raid/rcruces/data/tmp/sc_5k.npy')


# -----------------------------------------------------------------------------
# MEAN matrix SC from PNI 7T all subjects
# Calculate the group mean connectome
# sc_5k_mean = np.mean(SC_5k, axis=2)
sc_5k_mean = np.load('/host/yeatman/local_raid/rcruces/data/tmp/sc_5k_mean.npy')

# Calculate the Strenght of the Group Mean SC & log it for visualization purposes
A = np.log(np.sum(sc_5k_mean, axis=1))

# Replace nan and -inf
A[np.isnan(A)] = 0
A[np.isneginf(A)] = 0

# Plot the histogram
import matplotlib.pyplot as plt
plt.hist(A,bins=100)

# Plot the Group Mean SC Strenght
plot_hemispheres(i5_lh, i5_rh, array_name=A, cmap='rocket', nan_color=(0, 0, 0, 1),
  zoom=1.3, size=(900, 750), embed_nb=False, color_range=(7,9),
  color_bar='right')

# -----------------------------------------------------------------------------
# Script Diffusion Maps
def fslr5k_dm_lr(mtx, mask_5k, Ngrad=3, log=False, S=0, kernel=None):
    """
    Diffusion Map embedding from the SC or GD matrix.
    Calculates the LEFT diffusion map and then 
    Align RIGHT to LEFT with procrustis
    Use 
    
    log     :  (bool) False for GD gradients <<<< might remove this
    mask_5k :  (int) [0,1])
    """
    if log != True:
        mtx_log = mtx
    else:
        # log transform the connectome
        mtx_log = np.log(mtx)
    
    # Replace NaN with 0
    mtx_log[np.isnan(mtx_log)] = 0
    
    # Replace negative infinite with 0
    mtx_log[np.isneginf(mtx_log)] = 0
    
    # Replace infinite with 0
    mtx_log[~np.isfinite(mtx_log)] = 0
    
    # Replace 0 with epsilon
    mtx_log[mtx_log==0] = np.finfo(float).eps
    
    # Left and right mask
    indx_L = np.where(mask_5k[0:4842]==1)[0]
    indx_R = np.where(mask_5k[4842:9684]==1)[0]
    
    # Left and right SC
    mtx_L = mtx_log[0:4842, 0:4842]
    mtx_R = mtx_log[4842:9684, 4842:9684]
    
    # Slice the matrix
    mtx_L_masked = mtx_L[indx_L, :]
    mtx_L_masked = mtx_L_masked[:, indx_L]  
    mtx_R_masked = mtx_R[indx_R, :]
    mtx_R_masked = mtx_R_masked[:, indx_R]

    # mtx Left hemi
    mtx_L = GradientMaps(n_components=Ngrad, random_state=None, approach='dm', kernel=kernel)
    mtx_L.fit(mtx_L_masked, sparsity=S)

    # mtx Right hemi
    mtx_R = GradientMaps(n_components=Ngrad, alignment='procrustes', kernel=kernel); # align right hemi to left hemi
    mtx_R.fit(mtx_R_masked, sparsity=S, reference=mtx_L.gradients_)

    # Left and right gradients concatenated
    mtx_gradients = np.concatenate((mtx_L.gradients_, mtx_R.aligned_), axis=0)
    
    # Boolean mask
    mask_surf = mask_5k != 0
    
    # Get the index of the non medial wall regions
    indx = np.where(mask_5k==1)[0]

    # Map gradients to surface
    grad = [None] * Ngrad
    for i, g in enumerate(mtx_gradients.T[0:Ngrad,:]):
        # create a new array filled with NaN values
        g_nan = np.full(mask_surf.shape, np.nan)
        g_nan[indx] = g
        grad[i] = g_nan
    
    return(mtx_gradients, grad)

def fslr5k_dm(mtx, mask, Ngrad=3, Smooth=False, S=0.9, kernel='normalized_angle'):
    """Create the gradients from the MPC matrix"""
    # Cleanup before diffusion embeding
    mtx[~np.isfinite(mtx)] = 0
    mtx[np.isnan(mtx)] = 0
    mtx[mtx==0] = np.finfo(float).eps

    # Get the index of the non medial wall regions
    indx = np.where(mask==1)[0]

    # Slice the matrix
    mtx_masked = mtx[indx, :]
    mtx_masked = mtx_masked[:, indx]

    # Calculate the gradients
    gm = GradientMaps(n_components=Ngrad, random_state=None, approach='dm', kernel=kernel)
    gm.fit(mtx_masked, sparsity=S)
    
    # Map gradients to surface
    grad = [None] * Ngrad
    
    # Boolean mask
    mask_surf = mask != 0
    
    for i, g in enumerate(gm.gradients_.T[0:Ngrad,:]):
        
        # create a new array filled with NaN values
        g_nan = np.full(mask_surf.shape, np.nan)
        g_nan[indx] = g
        grad[i] = g_nan
    
    return(gm, grad)


# -----------------------------------------------------------------------------
# Diffusion maps and Plotting
Ngrad=10

# Calculate the gradients: FULL matrix
_, gradA = fslr5k_dm(sc_5k_mean, mask_new, Ngrad=Ngrad, S=0, kernel=None)

# plot the gradients
Nplot=5
labels=['EV'+str(x) for x in list(range(1,Nplot+1))]
plot_hemispheres(i5_lh, i5_rh, array_name=gradA[0:Nplot], cmap='RdBu_r', nan_color=(0, 0, 0, 1),
  zoom=1.3, size=(900, 750), embed_nb=False, color_range='sym',
  color_bar='right', label_text={'left': labels}, screenshot=False,
  filename='/home/bic/rcruces/Desktop/PsNI_SC-3G.png')

# Cleaner mask withoput outliers
mask_new = np.where(gradA[0]<0, 0, 1) * mask_5k


# Calculate the gradients: LEFT to RIGHT matrix
_, gradB = fslr5k_dm_lr(sc_5k_mean, mask_new, Ngrad=Ngrad, S=0, kernel=None, log=False)

# plot the gradients
plot_hemispheres(i5_lh, i5_rh, array_name=gradB[0:Nplot], cmap='RdBu_r', nan_color=(0, 0, 0, 1),
  zoom=1.3, size=(900, 750), embed_nb=False, color_range='sym',
  color_bar='right', label_text={'left': labels}, screenshot=False,
  filename='/home/bic/rcruces/Desktop/PNI_SC-3G.png')