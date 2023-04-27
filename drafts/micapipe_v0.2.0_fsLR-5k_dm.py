#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 08:32:03 2023

@author: rcruces
"""
# -------------------------------------------
# LIBRARIES
import os
import glob
import nibabel as nb
import numpy as np
from brainspace.mesh.mesh_io import read_surface
from brainspace.plotting import plot_hemispheres
from brainspace.mesh.array_operations import smooth_array
from brainspace.gradient import GradientMaps

# -------------------------------------------
# SET VARIABLES
micapipe='/data_/mica1/01_programs/micapipe-v0.2.0'
out = '/data_/mica3/BIDS_PNI/derivatives/micapipe_v0.2.0'
Ngrad = 3
dataset ='PNI'
Save_png='/home/bic/rcruces/Desktop/'
Save=False
# Set working directory
os.chdir(out)

# -------------------------------------------
# LOAD DATA
# fsLR-5k: inflated
c5k_lhi = read_surface(micapipe + '/surfaces/fsLR-5k.L.inflated.surf.gii', itype='gii')
c5k_rhi = read_surface(micapipe + '/surfaces/fsLR-5k.R.inflated.surf.gii', itype='gii')

# fsLR-5k: mask
mask_lh = nb.load(micapipe + '/surfaces/fsLR-5k.L.mask.shape.gii').darrays[0].data
mask_rh = nb.load(micapipe + '/surfaces/fsLR-5k.R.mask.shape.gii').darrays[0].data

# Concatenate both masks left and right
mask_fsRL10k = np.concatenate((mask_lh, mask_rh), axis=0)
# map_to_labels will use a boolean array as input
mask = mask_fsRL10k != 0

# Get the index of the non medial wall regions
indx = np.where(mask==1)[0]

# Number vertices per hemisphere
n5k = 4842

# -------------------------------------------
# DEFINE FUNCTIONS
def load_5k(File, Ndim, mod=None):
    # load the matrix
    mtx = np.loadtxt(File, dtype=float, delimiter=' ')
    
    # Mirror the matrix
    c5k = np.triu(mtx,1)+mtx.T
    
    return c5k

def load_5k_GD(File, Ndim, mod=None):
    # load the matrix
    mtx = np.loadtxt(File, dtype=float, delimiter=' ')
    return(mtx)

def smooth_surf(surf_l, surf_r, points, Kernel='uniform', Niter=20, Relax=0.5, Mask=mask):
    '''
    This function smooth an array (left and right) on a given surface
    Parameters
    ----------
    surf_l : np.array left surface 
    surf_r : np.array right surface
    points : np.array surface data
    Kernel : {'uniform', 'gaussian', 'inverse_distance'}
    Niter  : int, optional NUmber of smooth iterations
    Relax  : relax : float, optional relaxation facto
    Mask   : str or 1D ndarray, optional
    Returns
    -------
    sdata : smoothed numpy array
    
    '''
    n2dim = int(Mask.shape[0]/2)
    Mask_l=Mask[0:n2dim]
    Mask_r=Mask[n2dim:Mask.shape[0]]
    sdata = np.concatenate((smooth_array(surf_l, points[0:4842],kernel=Kernel, n_iter=Niter,relax=Relax, mask=Mask_l), 
                          smooth_array(surf_r, points[4842:4842*2],kernel=Kernel, n_iter=Niter,relax=Relax, mask=Mask_r)), axis=0)
    return(sdata)

# Load files
def load_connectomes(files, Ndim, func, Print=False):
    # Load all the matrices
    M=np.empty([Ndim*2, Ndim*2, len(files)], dtype=float)
    for i, f in enumerate(files):
        if Print==True: print(f)
        M[:,:,i] = func(f, Ndim)
    return M

# -------------------------------------------
# GD - Geodesic Distance

# Load connectomes
files=sorted(glob.glob('sub-*/ses-*/dist/*_surf-fsLR-5k_GD.txt'))
GD = load_connectomes(files, n5k, load_5k_GD, Print=True)

# --------------------------------------------------------------------------------------
# TEST save matrices as gii and the upload them again (check the loading times)
# txt takes forever :(
def save_gii(data_array, file_name):
    # Initialize gifti: NIFTI_INTENT_SHAPE - 2005, FLOAT32 - 16
    gifti_data = nb.gifti.GiftiDataArray(data=data_array, intent=2005, datatype=16)

    # this is the GiftiImage class
    gifti_img = nb.gifti.GiftiImage(meta=None, darrays=[gifti_data])

    # Save the new GIFTI file
    nb.save(img=gifti_img, filename=file_name)

# Save each matrix as gifti
for i,n, in enumerate(files):
    print(i+n)

# --------------------------------------------------------------------------------------

# Mean matrix across the z axis (subjects)
gd_mean = np.mean(GD, axis=2)

# left and right GD
GD_l=gd_mean[0:n5k, 0:n5k]
GD_r=gd_mean[n5k:n5k*2, n5k:n5k*2]

# Get the index of the non medial wall regions
indx_h = np.where(mask_lh==1)[0]

# Slice the matrix
GDl_m = GD_l[indx_h, :]
GDl_m = GDl_m[:, indx_h]

# Slice the matrix
GDr_m = GD_r[indx_h, :]
GDr_m = GDr_m[:, indx_h]

# Diffusion Embedding mapping
S=0.8
gm_GD = GradientMaps(n_components=Ngrad, random_state=None, approach='dm', kernel='normalized_angle')
gm_GD.fit(GDl_m, sparsity=S)
gm_GD_r = GradientMaps(n_components=Ngrad, alignment='procrustes', kernel='normalized_angle'); # align right hemi to left hemi
gm_GD_r.fit(GDr_m, sparsity=S, reference=gm_GD.gradients_)

# Left and right gradients concatenated
GD_gradients = np.concatenate((gm_GD.gradients_, gm_GD_r.aligned_), axis=0)

# Map gradients to surface
grad = [None] * Ngrad
for i, g in enumerate(GD_gradients.T[0:Ngrad,:]):
    g_nan = np.full(mask.shape, np.nan)
    g_nan[indx] = g
    # fill in the calculated values into the corresponding indices of the new array
    grad[i] = smooth_surf(c5k_lhi, c5k_rhi, g_nan, Niter=3, Relax=0.3)
    grad[i] = g_nan

# plot the gradients
labels=['G'+str(x) for x in list(range(1,Ngrad+1))]
plot_hemispheres(c5k_lhi, c5k_rhi, array_name=grad, cmap='RdBu_r', nan_color=(0, 0, 0, 1),
  zoom=1.3, size=(900, 750), embed_nb=True, color_range='sym', color_bar='right', label_text={'left': labels},
  screenshot=Save, filename=Save_png + dataset + '_GD_dm_fsLR-5k_MICs.png')  

# Plot the mean GD for QC
feat=np.sum(gd_mean,axis=1)
Range=(np.quantile(feat, 0.1), np.quantile(feat, 0.99))
plot_hemispheres(c5k_lhi, c5k_rhi, array_name=feat, cmap='Blues', nan_color=(0, 0, 0, 1),
                      zoom=1.3, size=(900, 750), embed_nb=True,
                      color_bar='right', layout_style='grid', color_range=Range,
                      label_text={'left': ['Lateral', 'Medial'], 'top': ['Left', 'Right']},
                      screenshot=Save, filename=Save_png + dataset + '_GD_fsLR-5k.png')

# -------------------------------------------
# MPC - Microstructure Profile Covariance
# Load all the connectomes
acq='acq-T1map'
files=sorted(glob.glob('sub-*/ses-*/mpc/'+acq+'/*_surf-fsLR-5k_desc-MPC.txt'))
CN = load_connectomes(files, 4842, load_5k)

# Mean matrix across the z axis (subjects)
cn_mean = np.mean(CN, axis=2)

# Cleanup before diffusion embeding
cn_mean[~np.isfinite(cn_mean)] = 0
cn_mean[np.isnan(cn_mean)] = 0
cn_mean[cn_mean==0] = np.finfo(float).eps

# Plot the mean MPC for QC
feat_s=smooth_surf(c5k_lhi, c5k_rhi,np.sum(cn_mean, axis=1), Niter=3, Relax=0.3)
Range=(np.quantile(feat_s, 0.01), np.quantile(feat_s, 0.95))
plot_hemispheres(c5k_lhi, c5k_rhi, array_name=feat_s, cmap='mako_r', nan_color=(0, 0, 0, 1),
                      zoom=1.3, size=(900, 750), embed_nb=True,
                      color_bar='right', layout_style='grid', color_range=Range,
                      label_text={'left': ['Lateral', 'Medial'], 'top': ['Left', 'Right']},
                      screenshot=False, filename=Save_png + dataset + '_MPC-' + acq + '_fsLR-5k.png')

# Slice the matrix
MPC_masked = cn_mean[indx, :]
MPC_masked = MPC_masked[:, indx]
MPC_masked.shape

# Calculate the gradients
MPCgm = GradientMaps(n_components=Ngrad, random_state=None, approach='dm', kernel='normalized_angle')
MPCgm.fit(MPC_masked, sparsity=0.9)

# Map gradients to surface
# other color 'RdYlBu_r'
grad = [None] * Ngrad
for i, g in enumerate(MPCgm.gradients_.T[0:Ngrad,:]):
    # create a new array filled with NaN values
    g_nan = np.full(mask.shape, np.nan)
    g_nan[indx] = g

    # fill in the calculated values into the corresponding indices of the new array
    grad[i] = smooth_surf(c5k_lhi, c5k_rhi, g_nan, Niter=3, Relax=0.3)

# plot the gradients
labels=['G'+str(x) for x in list(range(1,Ngrad+1))]
plot_hemispheres(c5k_lhi, c5k_rhi, array_name=grad, cmap='RdBu_r', nan_color=(0, 0, 0, 1),
  zoom=1.3, size=(900, 750), embed_nb=True, color_range='sym', color_bar='right', label_text={'left': labels},
  screenshot=False, filename=Save_png + dataset + '_MPC-' + acq + '_dm_fsLR-5k.png')  

# -------------------------------------------
# FC - Functional Connectome

# -------------------------------------------
# SC - Structural Connectome

# -------------------------------------------
# 