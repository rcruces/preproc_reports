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
import seaborn
import cmocean
from brainspace.mesh.mesh_io import read_surface
from brainspace.plotting import plot_hemispheres
from brainspace.mesh.array_operations import smooth_array
from brainspace.gradient import GradientMaps
cmaps = cmocean.cm.cmap_d

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
def load_5k(File, Ndim):
    # load the matrix
    #mtx = np.loadtxt(File, dtype=float, delimiter=' ')
    mtx = nb.load(File).darrays[0].data
    # Mirror the matrix
    c5k = np.triu(mtx,1)+mtx.T
    
    return c5k

def load_5k_GD(File, Ndim):
    # load the matrix
    #mtx = np.loadtxt(File, dtype=float, delimiter=' ')
    mtx = nb.load(File).darrays[0].data
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

# --------------------------------------------------------------------------------------
# GD - Geodesic Distance
# --------------------------------------------------------------------------------------

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
for i,f, in enumerate(files):
    print(f)
    save_gii(GD[:,:,i], f.replace('txt', 'shape.gii'))

files=sorted(glob.glob('sub-*/ses-*/dist/*_surf-fsLR-5k_GD.shape.gii'))
GD = load_connectomes(files, n5k, load_5k_GD, Print=True)
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
    grad[i] = g_nan

# plot the gradients
labels=['G'+str(x) for x in list(range(1,Ngrad+1))]
plot_hemispheres(c5k_lhi, c5k_rhi, array_name=grad, cmap='RdBu_r', nan_color=(0, 0, 0, 1),
  zoom=1.3, size=(900, 750), embed_nb=True, color_range='sym', color_bar='right', label_text={'left': labels},
  screenshot=Save, filename=Save_png + dataset + 'fsLR-5k_GD-dm.png')  

# Plot the mean GD for QC
feat=np.sum(gd_mean,axis=1)
Range=(np.quantile(feat, 0.1), np.quantile(feat, 0.99))
plot_hemispheres(c5k_lhi, c5k_rhi, array_name=feat, cmap='Blues', nan_color=(0, 0, 0, 1),
                      zoom=1.3, size=(900, 750), embed_nb=True,
                      color_bar='right', layout_style='grid', color_range=Range,
                      label_text={'left': ['Lateral', 'Medial'], 'top': ['Left', 'Right']},
                      screenshot=Save, filename=Save_png + dataset + 'fsLR-5k_GD-colsum.png')

# --------------------------------------------------------------------------------------
# MPC - Microstructure Profile Covariance
# --------------------------------------------------------------------------------------
def dm_5k(CN, mask, modality=None, Ngrad=3, S=0.9, Smooth=False,
          plot_group=True, mod='', qRange=(0.01,0.95),
          cmap='cmo.deep'):
    # Number of vertices per one hemisphere
    #N = 4842
    if modality=="GD":
        print("[Info]... GD DM will be calculated per hemisphere")
    elif modality=="SC":
        print("[Info]... SC DM will be calculated per hemisphere")

    # Mean matrix across the z axis (subjects)
    cn_mean = np.mean(CN, axis=2)
    
    # Cleanup before diffusion embeding
    cn_mean[~np.isfinite(cn_mean)] = 0
    cn_mean[np.isnan(cn_mean)] = 0
    cn_mean[cn_mean==0] = np.finfo(float).eps
    
    # Get the index of the non medial wall regions
    indx = np.where(mask==1)[0]

    # Slice the matrix
    CN_masked = cn_mean[indx, :]
    CN_masked = CN_masked[:, indx]

    print("[INFO]... Calculating the gradients")
    gm = GradientMaps(n_components=Ngrad, random_state=None, approach='dm', kernel='normalized_angle')
    gm.fit(CN_masked, sparsity=S)
    
    # Map gradients to surface
    grad = [None] * Ngrad
    for i, g in enumerate(gm.gradients_.T[0:Ngrad,:]):
        # create a new array filled with NaN values
        g_nan = np.full(mask.shape, np.nan)
        g_nan[indx] = g

        # fill in the calculated values into the corresponding indices of the new array
        if Smooth==True:
            grad[i] = smooth_surf(c5k_lhi, c5k_rhi, g_nan, Niter=3, Relax=0.3)
        else:
            grad[i] = g_nan
    
    # plot the gradients
    labels=['G'+str(x) for x in list(range(1,Ngrad+1))]
    p=plot_hemispheres(c5k_lhi, c5k_rhi, array_name=grad, cmap='RdBu_r', nan_color=(0, 0, 0, 1),
      zoom=1.3, size=(900, 750), embed_nb=False, color_range='sym', color_bar='right', label_text={'left': labels},
      screenshot=Save, filename=Save_png + dataset + 'fsLR-5k_'+mod+'_dm.png')  
    
    return(gm,p)

def dm_5k_plot(gm, mask, Smooth=True, plot_group=True,
               mod='', qRange=(0.01,0.95), cmap='cmo.deep'):
    # Map gradients to surface
    grad = [None] * Ngrad
    for i, g in enumerate(gm.gradients_.T[0:Ngrad,:]):
        # create a new array filled with NaN values
        g_nan = np.full(mask.shape, np.nan)
        g_nan[indx] = g

        # fill in the calculated values into the corresponding indices of the new array
        if Smooth==True:
            grad[i] = smooth_surf(c5k_lhi, c5k_rhi, g_nan, Niter=3, Relax=0.3)
        else:
            grad[i] = g_nan
    
    # plot the gradients
    labels=['G'+str(x) for x in list(range(1,Ngrad+1))]
    plot_hemispheres(c5k_lhi, c5k_rhi, array_name=grad, cmap='RdBu_r', nan_color=(0, 0, 0, 1),
      zoom=1.3, size=(900, 750), embed_nb=True, color_range='sym', color_bar='right', label_text={'left': labels},
      screenshot=Save, filename=Save_png + dataset + 'fsLR-5k_'+mod+'_dm.png')  
    
    if plot_group==True:
        # Plot the group mean Connectome
        if Smooth==True:
            feat=smooth_surf(c5k_lhi, c5k_rhi,np.sum(cn_mean, axis=1), Niter=3, Relax=0.3)
        else:
            feat=np.sum(cn_mean, axis=1)
        Range=(np.quantile(feat, qRange[0]), np.quantile(feat, qRange[0]))
        plot_hemispheres(c5k_lhi, c5k_rhi, array_name=feat*mask, cmap=cmap, nan_color=(0, 0, 0, 1),
                         zoom=1.3, size=(900, 750), embed_nb=True,
                         color_bar='right', layout_style='grid', color_range=Range,
                         label_text={'left': ['Lateral', 'Medial'], 'top': ['Left', 'Right']},
                         screenshot=Save, filename=Save_png + dataset + 'fsLR-5k_'+mod+'_colsum.png')

    

# Load all the connectomes
acq='acq-T1map'
files=sorted(glob.glob('sub-*/ses-*/mpc/'+acq+'/*_surf-fsLR-5k_desc-MPC.shape.gii'))
CN = load_connectomes(files, 4842, load_5k)

# Mean matrix across the z axis (subjects)
cn_mean = np.mean(CN, axis=2)

feat=np.abs(np.sum(cn_mean, axis=1))
Range=(np.quantile(feat, 0.001), np.quantile(feat, 0.1))
plot_hemispheres(c5k_lhi, c5k_rhi, array_name=feat, cmap='rocket', nan_color=(0, 0, 0, 1),
                 zoom=1.3, size=(900, 750), embed_nb=False,
                 color_bar='right', layout_style='grid', color_range=Range,
                 label_text={'left': ['Lateral', 'Medial'], 'top': ['Left', 'Right']},
                 screenshot=True, filename=Save_png + dataset + 'fsLR-5k_MPC-qT1_colsum_mask.png')

# Save each matrix as gifti
for i,f, in enumerate(files):
    print(f)
    save_gii(CN[:,:,i], f.replace('txt', 'shape.gii'))

dm_5k(CN, mask, modality=None, Ngrad=3, S=0.9, Smooth=False,
          plot_group=True, mod='', qRange=(0.01,0.95),
          cmap='cmo.deep')

# ------------------------------------------------------------
acq='acq-MTR'
files=sorted(glob.glob('sub-*/ses-*/mpc/'+acq+'/*_surf-fsLR-5k_desc-MPC.shape.gii'))
mtr = load_connectomes(files, 4842, load_5k, Print=True)

# Save each matrix as gifti
for i,f, in enumerate(files):
    print(f)
    save_gii(mtr[:,:,i], f.replace('txt', 'shape.gii'))

# Mean matrix across the z axis (subjects)
cn_mean = np.mean(mtr, axis=2)

# Cleanup before diffusion embeding
cn_mean[~np.isfinite(cn_mean)] = 0
cn_mean[np.isnan(cn_mean)] = 0
cn_mean[cn_mean==0] = np.finfo(float).eps

feat=smooth_surf(c5k_lhi, c5k_rhi,np.sum(cn_mean, axis=1), Niter=3, Relax=0.3)

feat=np.abs(np.sum(cn_mean, axis=1))
Range=(np.quantile(feat, 0.000001), np.quantile(feat, 0.2))
plot_hemispheres(c5k_lhi, c5k_rhi, array_name=feat, cmap='rocket', nan_color=(0, 0, 0, 1),
                 zoom=1.3, size=(900, 750), embed_nb=False,
                 color_bar='right', layout_style='grid', color_range=Range,
                 label_text={'left': ['Lateral', 'Medial'], 'top': ['Left', 'Right']},
                 screenshot=True, filename=Save_png + dataset + 'fsLR-5k_MPC-MTR_colsum.png')

# Get the index of the non medial wall regions
indx = np.where(mask==1)[0]

# Slice the matrix
CN_masked = cn_mean[indx, :]
CN_masked = CN_masked[:, indx]

dm_5k(mtr, mask, modality=None, Ngrad=3, S=0.9, Smooth=False,
          plot_group=True, mod='MPC-MTR', qRange=(0.01,0.95),
          cmap='cmo.deep')

# ------------------------------------------------------------
acq='acq-T2star'
files=sorted(glob.glob('sub-*/ses-*/mpc/'+acq+'/*_surf-fsLR-5k_desc-MPC.txt'))
t2s = load_connectomes(files, 4842, load_5k, Print=True)

# Save each matrix as gifti
for i,f, in enumerate(files):
    print(f)
    save_gii(t2s[:,:,i], f.replace('txt', 'shape.gii'))
# ------------------------------------------------------------

# Mean matrix across the z axis (subjects)
cn_mean = np.mean(t2s, axis=2)

# Cleanup before diffusion embeding
cn_mean[~np.isfinite(cn_mean)] = 0
cn_mean[np.isnan(cn_mean)] = 0
cn_mean[cn_mean==0] = np.finfo(float).eps

feat=smooth_surf(c5k_lhi, c5k_rhi,np.sum(cn_mean, axis=1), Niter=3, Relax=0.3)
Range=(np.quantile(feat, 0.1), np.quantile(feat, 0.95))
plot_hemispheres(c5k_lhi, c5k_rhi, array_name=feat*mask, cmap='rocket', nan_color=(0, 0, 0, 1),
                 zoom=1.3, size=(900, 750), embed_nb=False,
                 color_bar='right', layout_style='grid', color_range=Range,
                 label_text={'left': ['Lateral', 'Medial'], 'top': ['Left', 'Right']},
                 screenshot=Save, filename=Save_png + dataset + 'fsLR-5k_MPC-MTR_colsum.png')

# Get the index of the non medial wall regions
indx = np.where(mask==1)[0]

# Slice the matrix
CN_masked = cn_mean[indx, :]
CN_masked = CN_masked[:, indx]

dm_5k(t2s, mask, modality=None, Ngrad=3, S=0.9, Smooth=True,
          plot_group=True, mod='', qRange=(0.01,0.95),
          cmap='cmo.deep')

# -------------------------------------------
# FC - Functional Connectome

# -------------------------------------------
# SC - Structural Connectome

# -------------------------------------------
# 