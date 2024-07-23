#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
7T precision imaging

FIGURE 2

@author: rcruces
"""

# -----------------------------------------------------------------------------
# Libraries
import pandas as pd
import numpy as np
import os
import glob
import nibabel as nb
import scipy.stats
from brainspace.mesh.mesh_io import read_surface
from brainspace.plotting import plot_hemispheres
from brainspace.gradient import GradientMaps
from brainspace.mesh.array_operations import smooth_array
from brainstat.stats.SLM import SLM
from brainstat.stats.terms import FixedEffect
from brainstat.datasets.base import combine_surfaces
import matplotlib.cm as cm
import cmocean.cm as cmo
import seaborn as sns
from brainspace.datasets import load_mask
import seaborn as sns
import matplotlib.pyplot as plt
import cmocean
import scipy.stats as st
cmaps = cmocean.cm.cmap_d

# -----------------------------------------------------------------------------
# Variables and paths
# Path to MICAPIPE
micapipe=os.popen("echo $MICAPIPE").read()[:-1]

# Load native mid surface
inf_lh = read_surface(micapipe + '/surfaces/fsLR-5k.L.inflated.surf.gii', itype='gii')
inf_rh = read_surface(micapipe + '/surfaces/fsLR-5k.R.inflated.surf.gii', itype='gii')
mask_lh =  nb.load(micapipe + '/surfaces/fsLR-5k.R.mask.shape.gii').darrays[0].data
mask_rh =  nb.load(micapipe + '/surfaces/fsLR-5k.L.mask.shape.gii').darrays[0].data
mask_10k = np.concatenate((mask_lh, mask_rh), axis=0)

# Labels and boolean mask
mask_surf = mask_10k != 0

# load stuff
mask_32k = load_mask(join=True)
inf32_lh = read_surface(micapipe + '/surfaces/fsLR-32k.L.inflated.surf.gii', itype='gii')
inf32_rh = read_surface(micapipe + '/surfaces/fsLR-32k.R.inflated.surf.gii', itype='gii')

# Shape of the fsLR-5k matrices
N5k = 9684

# Set dataset PNI as working directory
os.chdir('/data_/mica3')

# List of subjects
subs_pni = ['PNC0'+str(f"{x:02}") for x in list(range(1,11))]

# full path to 
dir_pni = 'BIDS_PNI/data_release'

# For PNI subjects
dirs_pni = []
for sub in subs_pni:
    dirs_pni.extend(glob.glob(f'{dir_pni}/derivatives/micapipe_v0.2.0/sub-{sub}/ses-*'))

# Sort the directories
dirs_pni = sorted(dirs_pni)

# -----------------------------------------------------------------------------
# Functions
def plot_connectome(mtx, Title='matrix plot', xlab='X', ylab='Y', col='rocket', vmin=None, vmax=None,
                   xticklabels='auto', yticklabels='auto',xrot=90, yrot=0, size=(15,10)):
    
    '''
    This optional function, only plots a connectome as a heatmap
    Parameters
    ----------
    mtx : np.array
    Returns
    -------
    f : plot
    '''
    f, ax = plt.subplots(figsize=size)
    g = sns.heatmap(mtx, ax=ax, cmap=col, vmin=vmin, vmax=vmax, xticklabels=xticklabels, yticklabels=yticklabels)
    g.set_xlabel(xlab)
    g.set_ylabel(ylab)
    g.set_title(Title)
    # Rotate the x-axis labels
    # rotate tick labels
    ax.set_xticklabels(ax.get_xticklabels(), rotation=xrot, ha='right')
    ax.set_yticklabels(ax.get_yticklabels(), rotation=yrot, ha='right')
    
def load_fc(File):
    """Loads and process a functional connectome"""
    
    # load the matrix
    FC = nb.load(File).darrays[0].data
    
    # Fisher transform
    FCz = np.arctanh(FC)

    # replace inf with 0
    FCz[~np.isfinite(FCz)] = 0
    
    # Mirror the matrix
    FCz = np.triu(FCz,1)+FCz.T
    return FCz

def get_id(str_dir):
    # Split the string by '/'
    parts = str_dir.split('/')

    # Extract the subject and session parts
    subject_part = parts[4].split('-')[-1]  # 'PNC001'
    session_part = parts[5].split('-')[-1]  # '01'

    # Concatenate them to get the desired result
    id_str = subject_part + 's' + session_part
    
    # result
    return(id_str)

def upper_tri_indexing(A):
    m = A.shape[0]
    r,c = np.triu_indices(m,1)
    return A[r,c]

# Smooth intencities
def smooth_surf(surf_l, surf_r, points, Mask, Kernel='uniform', Niter=3, Relax=0.35):
    '''
    This function smooth an array (left and right) on a given surface
    Parameters
    ----------
    surf_l : np.array left surface 
    surf_r : np.array right surface
    points : np.array surface data
    Kernel : {'uniform', 'gaussian', 'inverse_distance'}
    Niter  : int, optional Number of smooth iterations
    Relax  : relax : float, optional relaxation factor
    Mask   : str or 1D ndarray, optional
    Returns
    -------
    sdata : smoothed numpy array
    
    '''
    Ndim = Mask.shape[0] 
    n2dim = int(Ndim/2)
    Mask_l=Mask[0:n2dim]
    Mask_r=Mask[n2dim:Ndim]
    sdata = np.concatenate((smooth_array(surf_l, points[0:n2dim],kernel=Kernel, n_iter=Niter,relax=Relax, mask=Mask_l), 
                          smooth_array(surf_r, points[n2dim:Ndim],kernel=Kernel, n_iter=Niter,relax=Relax, mask=Mask_r)), axis=0)
    return(sdata)

def hist_cmap(surf_array, lincmap=cm.bone, b=100, cmap_qt=[0.01,0.99],
              Title='Histogram', xlab='x-label', crange=None, x_range=None):
    
    # Remove NaN from the variance
    array_ft = surf_array[~np.isnan(surf_array)]
    n, bins = np.histogram(array_ft, bins=b)  # Compute the histogram data
    bin_centers = (bins[:-1] + bins[1:]) / 2  # Compute the bin centers
    if crange == None:
        colored_bins = lincmap(np.interp(bin_centers, [np.quantile(surf_array, cmap_qt[0] ), 
                                                       np.quantile(surf_array, cmap_qt[1] )], [0, 1]))
    else:
        colored_bins = lincmap(np.interp(bin_centers, [crange[0], crange[1] ], [0, 1]))
        
    # Create a figure with a transparent background
    fig, ax = plt.subplots()
    fig.patch.set_alpha(0.0)  # Set the background of the figure to be transparent
    ax.patch.set_alpha(0.0)  # Set the background of the axes to be transparent
    
    # Plot the histogram with the desired colormap
    plt.bar(bin_centers, n, width=np.diff(bins), color=colored_bins, alpha=1)
    
    # Remove upper and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Customize ticks and labels
    plt.xlabel(xlab, fontsize=10)  # Increase font size for x-label
    plt.ylabel('Frequency', fontsize=10)  # Increase font size for y-label
    plt.title(Title, fontsize=12)  # Increase font size for title
    
    # Set x-axis range if provided
    if x_range is not None:
        plt.xlim(x_range)
    
    plt.show()

# -----------------------------------------------------------------------------
# Group DM
# Create the gradients from a matrix
def array_dm(MPC, mpc_mask, Ngrad=3, Smooth=False, kernel='normalized_angle', reference=None, S=0.9, r2l=False):
    # Cleanup before diffusion embeding
    MPC[~np.isfinite(MPC)] = 0
    MPC[np.isnan(MPC)] = 0
    MPC[MPC==0] = np.finfo(float).eps

    # Get the index of the non medial wall regions
    indx = np.where(mpc_mask==1)[0]

    # Slice the matrix
    MPC_masked = MPC[indx, :]
    MPC_masked = MPC_masked[:, indx]
    mask_shape = MPC_masked.shape[1]
    mask_shape2 = int(mask_shape/2)

    # Calculate the gradients with no reference (procrustis alignment)
    if reference == None:
        # Align Right to Left
        if r2l == True:
            # Left hemi
            dgm_L = GradientMaps(n_components=Ngrad, random_state=None, approach='dm', kernel='normalized_angle')
            dgm_L.fit(MPC_masked[0:mask_shape2, 0:mask_shape2], sparsity=S)
            
            # Right hemi
            dgm_R = GradientMaps(n_components=Ngrad, alignment='procrustes', kernel='normalized_angle'); # align right hemi to left hemi
            dgm_R.fit(MPC_masked[mask_shape2:mask_shape, mask_shape2:mask_shape], sparsity=S, reference=dgm_L.gradients_)
            
            # Return boths DM objects in one list
            dgm =[dgm_L, dgm_R]
            
            # Concatenate both gradients
            dm_ev = np.concatenate((dgm_L.gradients_, dgm_R.aligned_), axis=0).T
            
        else:
            dgm = GradientMaps(n_components=Ngrad, random_state=None, approach='dm', kernel=kernel)
            dgm.fit(MPC_masked, sparsity=S)
            dm_ev = dgm.gradients_.T
    else:
        # Align Right to Left
        if r2l == True:
            # Left hemi
            dgm_L = GradientMaps(n_components=Ngrad, random_state=None, approach='dm', alignment='procrustes', kernel='normalized_angle')
            dgm_L.fit(MPC_masked[0:mask_shape2, 0:mask_shape2], sparsity=S, reference=reference[0].gradients_)
            
            # Right hemi
            dgm_R = GradientMaps(n_components=Ngrad, random_state=None, approach='dm', alignment='procrustes', kernel='normalized_angle'); # align right hemi to left hemi
            dgm_R.fit(MPC_masked[mask_shape2:mask_shape, mask_shape2:mask_shape], sparsity=S, reference=reference[1].aligned_)
            
            # Return boths DM objects in one list
            dgm =[dgm_L, dgm_R]
            
            # Concatenate both gradients
            dm_ev = np.concatenate((dgm_L.aligned_, dgm_R.aligned_), axis=0).T
            
        else:
            dgm = GradientMaps(n_components=Ngrad, random_state=None, approach='dm', alignment='procrustes', kernel='normalized_angle'); # align to mean Gradient
            dgm.fit(MPC_masked, sparsity=S, reference=reference.gradients_)
            dm_ev = dgm.aligned_.T
    
    # Map gradients to surface
    grad = [None] * Ngrad
    for i, g in enumerate(dm_ev[0:Ngrad,:]):
        # create a new array filled with NaN values
        g_nan = np.full(mask_surf.shape, np.nan)
        g_nan[indx] = g

        # fill in the calculated values into the corresponding indices of the new array
        if Smooth==True:
            grad[i] = smooth_surf(inf_lh, inf_rh, g_nan, mask_surf,Niter=3, Relax=0.35, Kernel='uniform')
        else:
            grad[i] = g_nan
    
    return(dgm, grad)

def extract_correlations(corr_matrix, num_subjects, num_sessions):
    assert corr_matrix.shape == (num_subjects * num_sessions, num_subjects * num_sessions), \
        "The correlation matrix shape does not match the number of subjects and sessions."
    
    pairs_correlations = []
    all_other_correlations = []

    for i in range(num_subjects):
        for j in range(num_sessions - 1):
            index1 = i * num_sessions + j
            index2 = i * num_sessions + j + 1
            pairs_correlations.append(corr_matrix[index1, index2])

    for i in range(num_subjects * num_sessions):
        for j in range(i):
            if (i % num_sessions == j % num_sessions and i // num_sessions == j // num_sessions + 1) or \
               (i % num_sessions == j % num_sessions + 1 and i // num_sessions == j // num_sessions):
                continue
            all_other_correlations.append(corr_matrix[i, j])

    return np.array(pairs_correlations), np.array(all_other_correlations)

# -----------------------------------------------------------------------------
# Structura connectome
# -----------------------------------------------------------------------------
def load_sc(File):
    """Loads and process a structura connectome"""
    
    # load the matrix
    mtx_sc = nb.load(File).darrays[0].data
    
    # Mirror the matrix
    mtx_sc = np.triu(mtx_sc,1)+mtx_sc.T
    
    return mtx_sc

# List the files
sc_file = []
for path in dirs_pni:
    sc_file.extend(glob.glob(f'{path}/dwi/connectomes/*_surf-fsLR-5k_desc-iFOD2-40M-SIFT2_full-connectome.shape.gii'))
N = len(sc_file)
print(f"Numer of FC: {N}") 

# Loads all the MPC fsLR-5k matrices
sc_5k=np.empty([N5k, N5k, len(sc_file)], dtype=float)
for i, f in enumerate(sc_file):
    sc_5k[:,:,i] = load_sc(f)

# Calculate the Group Mean connectome {5k x 5k}
sc_5k_mean = np.mean(sc_5k, axis=2)

# -----------------------------------------------------------------------------
# Mean
# Log transfor the matrix for visualization
sc_5k_mean_log = np.log(sc_5k_mean)
sc_5k_mean_log[np.isneginf(sc_5k_mean_log)] = 0

# Calculate the mean of each column {10k x 1}
sc_5k_mean_log = np.mean(sc_5k_mean_log, axis=1)

# plot the column mean of the mean connectome surface
crange=(np.quantile(sc_5k_mean_log[mask_surf], 0.05), np.quantile(sc_5k_mean_log[mask_surf], 0.975))
sc_5k_mean_log[mask_surf == False] = np.nan

plot_hemispheres(inf_lh, inf_rh, array_name=sc_5k_mean_log, cmap='bone', nan_color=(0, 0, 0, 1), 
                 transparent_bg=True, zoom=1.3, size=(900, 750), embed_nb=True, color_bar='right', 
                 layout_style='grid', color_range=crange,
                 label_text={'left': ['Lateral', 'Medial'], 'top': ['Left', 'Right']}, screenshot=False)

# Histogram of the mean
hist_cmap(sc_5k_mean_log[mask_surf], Title='SC Group Mean', xlab='Mean SC', lincmap=cm.bone)

# -----------------------------------------------------------------------------
# Variance

# Calculate the mean colum FC for all subjects
sc_5k_sub = np.mean(sc_5k, axis=1)

# Calculate the variance
map_var = np.var(sc_5k_sub, axis=1)
#map_var = np.abs(map_mean/np.std(surf_map, axis=0))

# Mask the medium wall
map_var[mask_surf == False] = np.nan

# plot the column mean of the mean connectome surface
crange=(np.quantile(map_var[mask_surf], 0.01), np.quantile(map_var[mask_surf], 0.95))

# Plot it on the surface
plot_hemispheres(inf_lh, inf_rh, array_name=map_var, cmap='cmo.deep_r', nan_color=(0, 0, 0, 1), 
                 transparent_bg=True, zoom=1.3, size=(900, 750), embed_nb=True, color_bar='right', 
                 layout_style='grid', color_range=crange,
                 label_text={'left': ['Lateral', 'Medial'], 'top': ['Left', 'Right']}, screenshot=False)

hist_cmap(map_var[mask_surf], Title='qT1 Group variance', xlab='Mean qT1', lincmap=cmo.deep_r,
          cmap_qt=[0.01,0.95])

# -----------------------------------------------------------------------------
# SC Group DM
Ngrad=3

# Calculate the GROUP gradients
dm_sc, grad = array_dm(sc_5k_mean, mask_surf, Ngrad=Ngrad, Smooth=False, S=0, r2l=True)

# plot the gradients
labels=['G'+str(x) for x in list(range(1,Ngrad+1))]
plot_hemispheres(inf_lh, inf_rh, array_name=grad, cmap='cmo.tarn', nan_color=(0.6, 0.6, 0.6, 1),
  zoom=1.3, size=(900, 150*Ngrad), embed_nb=True, color_range='sym',
  color_bar='right', label_text={'left': labels}, screenshot=True,
  filename=f'/home/bic/rcruces/Desktop/7t_fmri/dm{Ngrad}-sc_group.png')

# -----------------------------------------------------------------------------
# Calculate each Subject's DM aligned to the group
Nsub = sc_5k.shape[2]

# Create an empty array to store the DM
dm_sub = np.empty((Ngrad, N5k, Nsub))

# Iterate over each subject
for i in range(0, Nsub):
    
    # Subject ID
    sub_id = get_id(sc_file[i])
    png_name=f'/home/bic/rcruces/Desktop/7t_fmri/dm{Ngrad}-sc_{sub_id}.png'
    print(f'{i} {sub_id}')
    
    # Calculate the DM
    _, dm_i = array_dm(sc_5k[:,:,i], mask_surf, Ngrad=Ngrad, Smooth=False, reference=dm_sc, r2l=True)
    dm_sub[:,:,i] = np.asarray(dm_i)
    
    labels=['G'+str(x) for x in list(range(1,Ngrad+1))]
    plot_hemispheres(inf_lh, inf_rh, array_name=dm_i, cmap='cmo.tarn', nan_color=(0.6, 0.6, 0.6, 1),
      zoom=1.3, size=(900, 150*Ngrad), embed_nb=True, color_range='sym',
      color_bar='right', label_text={'left': labels}, screenshot=True,
      filename=png_name)

# -----------------------------------------------------------------------------
# Connectome based similarity

# Calculate the mean colum FC for all subjects
sc5k_sub = np.mean(sc_5k, axis=1)

# Correlation
corr = np.corrcoef(sc5k_sub[mask_surf,:].T)

# Similarity measurements
reliability, uniformity = extract_correlations(corr, 10, 2)
np.mean(uniformity)
np.mean(reliability)

# Array for the histogram
corr_sym = upper_tri_indexing(corr)

bids_ids = [get_id(dir) for dir in sc_file]

# Get the indices of the upper triangle, excluding the diagonal
upper_tri_indices = np.triu_indices_from(corr, k=1)

# Replace the upper triangle values with NaN
corr[upper_tri_indices] = np.nan

plot_connectome(corr, 'Subject similarity - SC', xlab=None, ylab=None, col='cmo.matter', vmin=0, vmax=1,
                yticklabels=bids_ids, xticklabels=bids_ids)

# Histogram
hist_cmap(corr_sym, lincmap=cmo.matter, b=50, cmap_qt=[0,1], x_range=[0,1],
              Title='Histogram', xlab='Connectome similarity')

# -----------------------------------------------------------------------------
# G1 Correlation between subjects
corr = np.corrcoef(dm_sub[0,mask_surf,:].T)

# Similarity measurements
reliability, uniformity = extract_correlations(corr, 10, 2)
np.mean(uniformity)
np.mean(reliability)

# Array for the histogram
corr_sym = upper_tri_indexing(corr)

bids_ids = [get_id(dir) for dir in sc_file]

# Get the indices of the upper triangle, excluding the diagonal
upper_tri_indices = np.triu_indices_from(corr, k=1)

# Replace the upper triangle values with NaN
corr[upper_tri_indices] = np.nan

plot_connectome(corr, 'Subject similarity - SC', xlab=None, ylab=None, col='cmo.matter', vmin=0, vmax=1,
                yticklabels=bids_ids, xticklabels=bids_ids)

# Histogram
hist_cmap(corr_sym, lincmap=cmo.matter, b=10, crange=[0,1], x_range=[0,1],
              Title='Histogram', xlab='G1 similarity')

np.save('/host/yeatman/local_raid/rcruces/git_here/preproc_reports/drafts/DM_SC.npy', dm_sub)
dm_sub = np.load('/host/yeatman/local_raid/rcruces/git_here/preproc_reports/drafts/DM_SC.npy')


# -----------------------------------------------------------------------------
# Functional MRI
# -----------------------------------------------------------------------------
# List the files
fc_file = []
for path in dirs_pni:
    fc_file.extend(glob.glob(f'{path}/func/desc-me_task-rest_bold/surf/*_surf-fsLR-5k_desc-FC.shape.gii'))
N = len(fc_file)
print(f"Numer of FC: {N}") 

# Loads all the MPC fsLR-5k matrices
fc_5k=np.empty([N5k, N5k, len(fc_file)], dtype=float)
for i, f in enumerate(fc_file):
    fc_5k[:,:,i] = load_fc(f)

# Calculate the Group Mean connectome {5k x 5k}
fc_5k_mean = np.mean(fc_5k, axis=2)

# -----------------------------------------------------------------------------
# Mean
# Calculate the mean of each column {10k x 1}
fc5k_mean = np.mean(fc_5k_mean, axis=1)

# plot the column mean of the mean connectome surface
crange=(np.quantile(fc5k_mean[mask_surf], 0.01), np.quantile(fc5k_mean[mask_surf], 0.99))
fc5k_mean[mask_surf == False] = np.nan

plot_hemispheres(inf_lh, inf_rh, array_name=fc5k_mean, cmap='bone', nan_color=(0, 0, 0, 1), 
                 transparent_bg=True, zoom=1.3, size=(900, 750), embed_nb=True, color_bar='right', 
                 layout_style='grid', color_range=crange,
                 label_text={'left': ['Lateral', 'Medial'], 'top': ['Left', 'Right']}, screenshot=False)


hist_cmap(fc5k_mean[mask_surf], Title='FC Group Mean', xlab='Mean FC', lincmap=cm.bone)

# Calculate the GROUP gradients
Ngrad=3
dm_fc, grad = array_dm(fc_5k_mean, mask_surf, Ngrad=Ngrad, Smooth=False)

# plot the gradients
labels=['G'+str(x) for x in list(range(1,Ngrad+1))]
plot_hemispheres(inf_lh, inf_rh, array_name=grad, cmap='cmo.tarn', nan_color=(0.6, 0.6, 0.6, 1),
  zoom=1.3, size=(900, 150*Ngrad), embed_nb=True, color_range='sym',
  color_bar='right', label_text={'left': labels}, screenshot=True,
  filename=f'/home/bic/rcruces/Desktop/7t_fmri/dm{Ngrad}-fc_group.png')

# -----------------------------------------------------------------------------
# Calculate each Subject's DM aligned to the group
Nsub = fc_5k.shape[2]

# Create an empty array to store the DM
dm_sub = np.empty((Ngrad, N5k, Nsub))

# Iterate over each subject
for i in range(0, Nsub):
    
    # Subject ID
    sub_id = get_id(dirs_pni[i])
    png_name=f'/home/bic/rcruces/Desktop/7t_fmri/dm{Ngrad}-fc_{sub_id}.png'
    print(f'{i} {sub_id}')
    
    # Calculate the DM
    _, dm_i = array_dm(fc_5k[:,:,i], mask_surf, Ngrad=Ngrad, Smooth=False, reference=dm_fc)
    dm_sub[:,:,i] = np.asarray(dm_i)
    
    labels=['G'+str(x) for x in list(range(1,Ngrad+1))]
    plot_hemispheres(inf_lh, inf_rh, array_name=dm_i, cmap='cmo.tarn', nan_color=(0.6, 0.6, 0.6, 1),
      zoom=1.3, size=(900, 150*Ngrad), embed_nb=True, color_range='sym',
      color_bar='right', label_text={'left': labels}, screenshot=True,
      filename=png_name)

# -----------------------------------------------------------------------------
# Connectome based similarity

# Calculate the mean colum FC for all subjects
fc5k_sub = np.mean(fc_5k, axis=1)

# Correlation
corr = np.corrcoef(fc5k_sub[mask_surf,:].T)

# Similarity measurements
reliability, uniformity = extract_correlations(corr, 10, 3)
np.mean(uniformity)
np.mean(reliability)

# Array for the histogram
corr_sym = upper_tri_indexing(corr)

bids_ids = [get_id(dir) for dir in dirs_pni]

# Get the indices of the upper triangle, excluding the diagonal
upper_tri_indices = np.triu_indices_from(corr, k=1)

# Replace the upper triangle values with NaN
corr[upper_tri_indices] = np.nan

plot_connectome(corr, 'Subject similarity - FC', xlab=None, ylab=None, col='cmo.matter', vmin=0, vmax=1,
                yticklabels=bids_ids, xticklabels=bids_ids)

# Histogram
hist_cmap(corr_sym, lincmap=cmo.matter, b=50, cmap_qt=[0,1], x_range=[0,1],
              Title='Histogram', xlab='Connectome similarity')

# -----------------------------------------------------------------------------
dm_sub = np.load('/host/yeatman/local_raid/rcruces/git_here/preproc_reports/drafts/DM_FC.npy')

# G1 Correlation between subjects
corr = np.corrcoef(dm_sub[0,mask_surf,:].T)

# Similarity measurements
reliability, uniformity = extract_correlations(corr, 10, 3)
np.mean(uniformity)
np.mean(reliability)

# Array for the histogram
corr_sym = upper_tri_indexing(corr)

bids_ids = [get_id(dir) for dir in dirs_pni]

# Get the indices of the upper triangle, excluding the diagonal
upper_tri_indices = np.triu_indices_from(corr, k=1)

# Replace the upper triangle values with NaN
corr[upper_tri_indices] = np.nan

plot_connectome(corr, 'Subject similarity - FC', xlab=None, ylab=None, col='cmo.matter', vmin=0, vmax=1,
                yticklabels=bids_ids, xticklabels=bids_ids)

# Histogram
hist_cmap(corr_sym, lincmap=cmo.matter, b=50, crange=[0,1], x_range=[0,1],
              Title='Histogram', xlab='G1 similarity')

#np.save('/host/yeatman/local_raid/rcruces/git_here/preproc_reports/drafts/DM_FC.npy', dm_sub)

# -----------------------------------------------------------------------------
# D2 based similarity

# -----------------------------------------------------------------------------
# qT1 MPC
# -----------------------------------------------------------------------------
def load_mpc(File):
    """Loads and process a MPC"""

    # Load file
    mpc = nb.load(File).darrays[0].data
    
    # Mirror the lower triangle
    mpc = np.triu(mpc,1)+mpc.T
    
    # Replace infinite values with epsilon
    mpc[~np.isfinite(mpc)] = np.finfo(float).eps
    
    # Replace 0 with epsilon
    mpc[mpc==0] = np.finfo(float).eps
    
    # retun the MPC
    return(mpc)

# MPC T1map
qmri='T1map'
mpc_file = []
for path in dirs_pni:
    mpc_file.extend(glob.glob(f'{path}/mpc/acq-{qmri}/*surf-fsLR-5k_desc-MPC.shape.gii'))
N = len(mpc_file) 
print(f"Number of MPC: {N}") 

# Loads all the MPC fsLR-5k matrices
mpc_5k=np.empty([N5k, N5k, len(mpc_file)], dtype=float)
for i, f in enumerate(mpc_file):
    mpc_5k[:,:,i] = load_mpc(f)

# Mean group MPC colums
mpc_5k_mean = np.mean(mpc_5k, axis=2)

# Calculate the GROUP gradients
Ngrad=3
dm_t1, grad = array_dm(mpc_5k_mean, mask_surf, Ngrad=Ngrad, Smooth=True)

# plot the gradients
labels=['G'+str(x) for x in list(range(1,Ngrad+1))]
plot_hemispheres(inf_lh, inf_rh, array_name=grad, cmap='cmo.tarn', nan_color=(0.6, 0.6, 0.6, 1),
  zoom=1.3, size=(900, 150*Ngrad), embed_nb=True, color_range='sym',
  color_bar='right', label_text={'left': labels}, screenshot=False,
  filename=f'/home/bic/rcruces/Desktop/7t_fmri/dm{Ngrad}-t1_group.png')

# -----------------------------------------------------------------------------
# Calculate each Subject's DM aligned to the group
Nsub = mpc_5k.shape[2]

# Create an empty array to store the DM
t1_sub = np.empty((Ngrad, N5k, Nsub))

# Iterate over each subject
for i in range(0, Nsub):
    
    # Subject ID
    sub_id = get_id(dirs_pni[i])
    png_name=f'/home/bic/rcruces/Desktop/7t_fmri/dm{Ngrad}-t1_{sub_id}.png'
    print(f'{i} {sub_id}')
    
    # Calculate the DM
    _, dm_i = array_dm(mpc_5k[:,:,i], mask_surf, Ngrad=Ngrad, Smooth=True, reference=dm_t1)
    t1_sub[:,:,i] = np.asarray(dm_i)
    
    labels=['G'+str(x) for x in list(range(1,Ngrad+1))]
    plot_hemispheres(inf_lh, inf_rh, array_name=dm_i, cmap='cmo.tarn', nan_color=(0.6, 0.6, 0.6, 1),
      zoom=1.3, size=(900, 150*Ngrad), embed_nb=True, color_range='sym',
      color_bar='right', label_text={'left': labels}, screenshot=True,
      filename=png_name)

#np.save('/host/yeatman/local_raid/rcruces/git_here/preproc_reports/drafts/DM_T1.npy', t1_sub)
t1_sub = np.load('/host/yeatman/local_raid/rcruces/git_here/preproc_reports/drafts/DM_T1.npy')

# -----------------------------------------------------------------------------
# Mean
# Calculate the mean of each column {10k x 1}
fc5k_mean = np.mean(mpc_5k_mean, axis=1)

# plot the column mean of the mean connectome surface
crange=(np.quantile(fc5k_mean[mask_surf], 0.01), np.quantile(fc5k_mean[mask_surf], 0.999))
fc5k_mean[mask_surf == False] = np.nan

plot_hemispheres(inf_lh, inf_rh, array_name=fc5k_mean, cmap='bone', nan_color=(0, 0, 0, 1), 
                 transparent_bg=True, zoom=1.3, size=(900, 750), embed_nb=True, color_bar='right', 
                 layout_style='grid', color_range=crange,
                 label_text={'left': ['Lateral', 'Medial'], 'top': ['Left', 'Right']}, screenshot=False)

hist_cmap(fc5k_mean[mask_surf], Title='qT1 Group Mean', xlab='Mean FC', lincmap=cm.bone,
          cmap_qt=[0.01,0.999])

# -----------------------------------------------------------------------------
# Variance
# Calculate the mean colum FC for all subjects
t15k_sub = np.mean(mpc_5k, axis=1)

# Calculate the variance
map_var = np.var(t15k_sub, axis=1)
#map_var = np.abs(map_mean/np.std(surf_map, axis=0))

# Mask the medium wall
map_var[mask_surf == False] = np.nan

# plot the column mean of the mean connectome surface
crange=(np.quantile(map_var[mask_surf], 0.01), np.quantile(map_var[mask_surf], 0.995))

# Plot it on the surface
plot_hemispheres(inf_lh, inf_rh, array_name=map_var, cmap='cmo.deep_r', nan_color=(0, 0, 0, 1), 
                 transparent_bg=True, zoom=1.3, size=(900, 750), embed_nb=True, color_bar='right', 
                 layout_style='grid', color_range=crange,
                 label_text={'left': ['Lateral', 'Medial'], 'top': ['Left', 'Right']}, screenshot=False)

hist_cmap(map_var[mask_surf], Title='qT1 Group variance', xlab='Mean qT1', lincmap=cmo.deep_r,
          cmap_qt=[0.01,0.995])

# -----------------------------------------------------------------------------
# G1 Correlation between subjects
corr = np.corrcoef(t1_sub[0,mask_surf,:].T)

# Array for the histogram
corr_sym = upper_tri_indexing(corr)

bids_ids = [get_id(dir) for dir in dirs_pni]

# Get the indices of the upper triangle, excluding the diagonal
upper_tri_indices = np.triu_indices_from(corr, k=1)

# Replace the upper triangle values with NaN
corr[upper_tri_indices] = np.nan

plot_connectome(corr, 'Subject similarity - qT1', xlab=None, ylab=None, col='cmo.matter', vmin=0, vmax=1,
                yticklabels=bids_ids, xticklabels=bids_ids)

# Histogram
hist_cmap(corr_sym, lincmap=cmo.matter, b=50, crange=[0,1], x_range=[0,1],
              Title='Histogram', xlab='G1 similarity')

# -----------------------------------------------------------------------------
# Group gradient
# Calculate the mean colum FC for all subjects
fc5k_sub = np.mean(fc_5k, axis=1)

# Calculate the mean colum FC for all subjects
t15k_sub = np.mean(mpc_5k, axis=1)

# COncatenate the data excluding the medial wall
sub_cns = np.concatenate((fc5k_sub[mask_surf,:], t15k_sub[mask_surf,:]), axis=0)

# Subjects correlation
corr = np.corrcoef(sub_cns.T)

# Array for the histogram
corr_sym = upper_tri_indexing(corr)

bids_ids = [get_id(dir) for dir in dirs_pni]

# Get the indices of the upper triangle, excluding the diagonal
upper_tri_indices = np.triu_indices_from(corr, k=1)

# Replace the upper triangle values with NaN
corr[upper_tri_indices] = np.nan

plot_connectome(corr, 'Subject similarity - qT1', xlab=None, ylab=None, col='cmo.matter', vmin=0, vmax=1,
                yticklabels=bids_ids, xticklabels=bids_ids)

# Histogram
hist_cmap(corr_sym, lincmap=cmo.matter, b=50, crange=[0,1], x_range=[0,1],
              Title='Histogram', xlab='G1 similarity')

gm = GradientMaps(n_components=Ngrad, random_state=None, approach='dm', kernel='normalized_angle')
gm.fit(corr, sparsity=0.4)

# Plot the gradients
g1 = gm.gradients_[:, 0]
g2 = gm.gradients_[:, 1]
g3 = gm.gradients_[:, 2]

# Creating figure
fig = plt.figure(figsize=(7, 5))
ax = fig.add_subplot(111, projection="3d")

# Create a colormap from tab10
colormap = plt.get_cmap('tab10')

# Assign colors based on index, repeating every 10 colors, every three points with the same color
colors = [colormap((i // 3) % 10) for i in range(len(g1))]

# Creating plot
ax.scatter3D(g1, g2, g3, color=colors, s=70)
plt.title("Subjects space")

ax.set_xlabel('Grad 1')
ax.set_ylabel('Grad 2')
ax.set_zlabel('Grad 3')

# Remove the outer box lines
ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False

# Show plot
plt.show()