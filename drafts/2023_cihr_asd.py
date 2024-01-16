#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 21:51:27 2023

@author: rcruces
"""
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
import seaborn as sns

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

# Smooth each intencity
def smooth_intensities(int_profile):
    smoothed_i = np.copy(int_profile)
   
    # smooth each intensity
    for i in range(int_profile.shape[0]):
        smoothed_i[i,:] = smooth_surf(inf_lh, inf_rh, int_profile[i,:], mask_surf,Niter=5, Relax=0.5, Kernel='uniform')
   
    return(smoothed_i)


# Function to build the MPC from an intencity profile
def build_mpc(data, mask):
    # If no parcellation is provided, MPC will be computed vertexwise
    I = data

    # Calculate mean across columns, excluding mask and any excluded labels input
    I_M = np.nanmean(np.float32(np.where(mask, I, np.nan)), axis=1)

    # Get residuals of all columns (controlling for mean)
    I_resid = np.zeros(I.shape)
    for c in range(I.shape[1]):
        y = I[:,c]
        x = I_M
        slope, intercept, _, _, _ = scipy.stats.linregress(x,y)
        y_pred = intercept + slope*x
        I_resid[:,c] = y - y_pred

    # Calculate correlation coefficient of the intesities with residuals
    R = np.corrcoef(I_resid, rowvar=False)

    # Log transform
    MPC = 0.5 * np.log( np.divide(1 + R, 1 - R) )
    MPC[np.isnan(MPC)] = 0
    MPC[np.isinf(MPC)] = 0

    # CLEANUP: correct diagonal and round values to reduce file size
    # Replace all values in diagonal by zeros to account for floating point error
    for i in range(0,MPC.shape[0]):
            MPC[i,i] = 0
    
    # Output MPC, microstructural profiles, and problem nodes
    return (MPC)

# Create the gradients from the MPC matrix
def mpc_dm(MPC, mpc_mask, Ngrad=3, Smooth=False):
    # Cleanup before diffusion embeding
    MPC[~np.isfinite(MPC)] = 0
    MPC[np.isnan(MPC)] = 0
    MPC[MPC==0] = np.finfo(float).eps

    # Get the index of the non medial wall regions
    indx = np.where(mpc_mask==1)[0]

    # Slice the matrix
    MPC_masked = MPC[indx, :]
    MPC_masked = MPC_masked[:, indx]
    MPC_masked.shape

    # Calculate the gradients
    MPCgm = GradientMaps(n_components=Ngrad, random_state=None, approach='dm', kernel='normalized_angle')
    MPCgm.fit(MPC_masked, sparsity=0.9)
    
    # Map gradients to surface
    grad = [None] * Ngrad
    for i, g in enumerate(MPCgm.gradients_.T[0:Ngrad,:]):
        # create a new array filled with NaN values
        g_nan = np.full(mask_surf.shape, np.nan)
        g_nan[indx] = g

        # fill in the calculated values into the corresponding indices of the new array
        if Smooth==True:
            grad[i] = smooth_surf(inf_lh, inf_rh, g_nan, mask_surf,Niter=3, Relax=0.35, Kernel='uniform')
        else:
            grad[i] = g_nan
    
    return(MPCgm, grad)

def map_to_labels5k(mpc_sliced, mask):
    # Get the index of the non medial wall regions
    mask_indx = np.where(mask==1)[0]
    # map to the labels
    labels_5k = np.full(mask.shape, np.nan)
    labels_5k[mask_indx] = mpc_sliced
    return(labels_5k)

    
# Load files
def make_mpc5k(files):
    # Load all the matrices
    Ndim=9684
    M=np.empty([Ndim, Ndim, len(files)], dtype=float)
    for i, f in enumerate(files):
        print(f)
        int_profile = nb.load(f).darrays[0].data
        M[:,:,i] = build_mpc(smooth_intensities(int_profile), mask_surf)

    return M

# -----------------------------------------------------------------------------
# Set dataset PNI as working directory
os.chdir('/data_/mica3/BIDS_PNI/derivatives/micapipe_v0.2.0')

# Load thickness to make a better mask
lh_files = sorted(glob.glob('sub-PNA*/ses-01/maps/*_hemi-L_surf-fsLR-5k_label-thickness.func.gii'))
rh_files = sorted(glob.glob('sub-PNA*/ses-01/maps/*_hemi-R_surf-fsLR-5k_label-thickness.func.gii'))

# Load all the thickness data
Nth=np.concatenate((nb.load(lh_files[0]).darrays[0].data, nb.load(rh_files[0]).darrays[0].data), axis=0).shape[0]

surf_map=np.empty([len(lh_files), Nth], dtype=float)
for i, _ in enumerate(lh_files):
    surf_map[i,:] = np.hstack(np.concatenate((nb.load(lh_files[i]).darrays[0].data, nb.load(rh_files[i]).darrays[0].data), axis=0))
    # Mean matrix across the x axis (vertices)
    map_mean = np.mean(surf_map, axis=0)

mask_surf = np.where(map_mean<0.3,False, True)
plot_hemispheres(inf_lh, inf_rh, array_name=mask_surf, cmap='rocket_r', nan_color=(0, 0, 0, 1),
                              zoom=1.3, size=(900, 750), embed_nb=False,
                              color_bar='right', layout_style='grid', color_range=(0.5,1),
                              label_text={'left': ['Lateral', 'Medial'], 'top': ['Left', 'Right']},
                              screenshot=False)

# -----------------------------------------------------------------------------
# All subjects 
# MPC files of Authism spectrum Disorder
files_all = sorted(glob.glob('sub-*/ses-01/mpc/acq-T1map/*_surf-fsLR-5k_desc-intensity_profiles.shape.gii'))

# make the MPC
mpc_all = make_mpc5k(files_all)

# Mean matrix across the z axis (subjects)
mpc_all_avg = np.mean(mpc_all, axis=2)

# Calculate the gradients
Ngrad=3
MPC_all, grad = mpc_dm(mpc_all_avg, mask_surf, Ngrad=Ngrad, Smooth=False)
Gall = grad[0]

# plot the gradients
labels=['G'+str(x) for x in list(range(1,Ngrad+1))]
plot_hemispheres(inf_lh, inf_rh, array_name=grad, cmap='RdBu_r', nan_color=(0, 0, 0, 1),
  zoom=1.3, size=(900, 750), embed_nb=False, color_range='sym',
  color_bar='right', label_text={'left': labels}, screenshot=False,
  filename='/home/bic/rcruces/Desktop/MPC_ASD_3G.png')

# -----------------------------------------------------------------------------
# ASD Group G1 
plot_hemispheres(inf_lh, inf_rh, array_name=Gall, cmap='RdBu_r',
                  embed_nb=False,  label_text={'left':['G1']}, color_bar='right',
                  layout_style='grid', color_range=(-0.2, 0.2),
                  zoom=1.3, size=(900, 750), nan_color=(0, 0, 0, 1), 
                  screenshot=True, filename='/home/bic/rcruces/Desktop/MPC_ALL_G1.png')


# -----------------------------------------------------------------------------
# Authism Spectrum Disorder

# MPC files of Authism spectrum Disorder
files_asd = sorted(glob.glob('sub-PNA*/ses-01/mpc/acq-T1map/*_surf-fsLR-5k_desc-intensity_profiles.shape.gii'))

# make the MPC
mpc_asd = make_mpc5k(files_asd)

# get the subject ID
ids_asd = [elem.split('/')[0] + '_' + elem.split('/')[1] for elem in files_asd]

# Mean matrix across the z axis (subjects)
mpc_asd_avg = np.mean(mpc_asd, axis=2)

# Calculate the gradients
Ngrad=3
MPC_dm_asd, grad = mpc_dm(mpc_asd_avg, mask_surf, Ngrad=Ngrad, Smooth=False)
Gasd = grad[0]

# plot the gradients
labels=['G'+str(x) for x in list(range(1,Ngrad+1))]
plot_hemispheres(inf_lh, inf_rh, array_name=grad, cmap='Spectral_r', nan_color=(0, 0, 0, 1),
  zoom=1.3, size=(900, 750), embed_nb=False, color_range='sym',
  color_bar='right', label_text={'left': labels}, screenshot=False,
  filename='/home/bic/rcruces/Desktop/MPC_ASD_3G.png')

# -----------------------------------------------------------------------------
# ASD Group G1 
plot_hemispheres(inf_lh, inf_rh, array_name=Gasd, cmap='Spectral_r',
                  embed_nb=False,  label_text={'left':['G1']}, color_bar='right',
                  layout_style='grid', color_range=(-0.2, 0.2),
                  zoom=1.3, size=(900, 750), nan_color=(0, 0, 0, 1), 
                  screenshot=True, filename='/home/bic/rcruces/Desktop/MPC_ASD_G1.png')

# -----------------------------------------------------------------------------
# Healthy SUbjects
# MPC files of Healthy controls
files_hc = sorted(glob.glob('sub-PNC*/ses-01/mpc/acq-T1map/*_surf-fsLR-5k_desc-intensity_profiles.shape.gii'))
# exclude PNC004 and PNC011
files_hc = [item for item in files_hc if "PNC004" not in item and "PNC011" not in item]

# get the subject ID
ids_hc = [elem.split('/')[0] + '_' + elem.split('/')[1] for elem in files_hc]

# Make MPC per subject
mpc_hc = make_mpc5k(files_hc)

# Mean matrix across the z axis (subjects)
mpc_hc_avg = np.mean(mpc_hc, axis=2)

# Calculate the gradients
Ngrad=3
MPC_dm, grad = mpc_dm(mpc_hc_avg, mask_surf, Ngrad=Ngrad, Smooth=False)
Ghc = grad[0]

# plot the gradients
labels=['G'+str(x) for x in list(range(1,Ngrad+1))]
plot_hemispheres(inf_lh, inf_rh, array_name=grad, cmap='Spectral_r', nan_color=(0, 0, 0, 1),
  zoom=1.3, size=(900, 750), embed_nb=False, color_range='sym',
  color_bar='right', label_text={'left': labels}, screenshot=False,
  filename='/home/bic/rcruces/Desktop/MPC_HC_3G.png')


# -----------------------------------------------------------------------------
# Healthy Control Group G1 
plot_hemispheres(inf_lh, inf_rh, array_name=Ghc, cmap='Spectral_r',
                  embed_nb=False,  label_text={'left':['G1']}, color_bar='right',
                  layout_style='grid', color_range=(-0.2, 0.2),
                  zoom=1.3, size=(900, 750), nan_color=(0, 0, 0, 1), 
  screenshot=True, filename='/home/bic/rcruces/Desktop/MPC_HC_G1.png')


# -----------------------------------------------------------------------------

def mpc_cleanup(mpc_mtx, mask):
    # Cleanup before diffusion embeding
    mpc_mtx[~np.isfinite(mpc_mtx)] = 0
    mpc_mtx[np.isnan(mpc_mtx)] = 0
    mpc_mtx[mpc_mtx==0] = np.finfo(float).eps
    
    # Get the index of the non medial wall regions
    indx = np.where(mask==1)[0]
    
    # Slice the matrix
    MPC_masked = mpc_mtx[indx, :]
    MPC_masked = MPC_masked[:, indx]
    
    return(MPC_masked)

# Cohen's D
# Concatenate ASD and HC MPCs
mpc_all = np.concatenate((mpc_hc, mpc_asd), axis=2)
ids = ids_hc + ids_asd

# Index of the controls
indx =np.array(["PNC" in f for f in ids])

# Empty array for the aligned gradients
# ROI x Ngradients x Subjects
MPCdm = np.empty((9046, Ngrad, mpc_all.shape[2]), dtype=float)

# Calculate and align a gradient per subject
for k, subMPC in enumerate(ids):
    print(subMPC)
    # Create gradient model
    MPCsub = GradientMaps(n_components=Ngrad, alignment='procrustes', kernel='normalized_angle')
    # Align to HC mean
    MPCsub.fit(mpc_cleanup(mpc_all[:,:,k], mask_surf), sparsity=0.9, reference=MPC_all.gradients_)
    MPCdm[:,:,k] = MPCsub.aligned_

from numpy import std, mean, sqrt

#correct if the population S.D. is expected to be equal for the two groups.
def cohen_d(x,y):
    nx = len(x)
    ny = len(y)
    dof = nx + ny - 2
    return (mean(x) - mean(y)) / sqrt(((nx-1)*std(x, ddof=1) ** 2 + (ny-1)*std(y, ddof=1) ** 2) / dof)

# -----------------------------------------------------------------------------
# Gradient 1 cohen's D
[ids[x] for x in [1,6,12]]
[ids[x] for x in [18,19,20]]
indx_C=[1,6,12]
indx_A=[18,19,20]

#G1_cd = np.array([cohen_d(MPCdm[x,0,indx_A], MPCdm[x,0,indx_C]) for x in range(0,MPCdm[:,0,0].shape[0]) ])
#G2_cd = np.array([cohen_d(MPCdm[x,1,indx_A], MPCdm[x,1,indx_C]) for x in range(0,MPCdm[:,0,0].shape[0]) ])
#G3_cd = np.array([cohen_d(MPCdm[x,2,indx_A], MPCdm[x,2,indx_C]) for x in range(0,MPCdm[:,0,0].shape[0]) ])

# calculta Cohen's D
G1_cd = np.array([cohen_d(MPCdm[x,0,~indx], MPCdm[x,0,indx]) for x in range(0,MPCdm[:,0,0].shape[0]) ])
G2_cd = np.array([cohen_d(MPCdm[x,1,~indx], MPCdm[x,1,indx]) for x in range(0,MPCdm[:,0,0].shape[0]) ])
G3_cd = np.array([cohen_d(MPCdm[x,2,~indx], MPCdm[x,2,indx]) for x in range(0,MPCdm[:,0,0].shape[0]) ])

#Gcd = [map_to_labels5k(G1_cd, mask_surf), 
#       map_to_labels5k(G2_cd, mask_surf), 
#       map_to_labels5k(G3_cd, mask_surf)]

# Smooth the gradients
MPCdmSgrad = np.empty((9684, MPCdm.shape[1], MPCdm.shape[2]), float)
MPCdmS = np.copy(MPCdm)
for G in range(MPCdm.shape[1]):
    for sub in range(MPCdm.shape[2]):
        MPCdmSgrad[:,G,sub] = smooth_surf(inf_lh, inf_rh, map_to_labels5k(MPCdm[:,G,sub], mask_surf), mask_surf,Niter=5, Relax=0.7)

# calculta Cohen's D
G1_cd = np.array([cohen_d(MPCdmSgrad[x,0,~indx], MPCdmSgrad[x,0,indx]) for x in range(0,MPCdmSgrad[:,0,0].shape[0]) ])
G2_cd = np.array([cohen_d(MPCdmSgrad[x,1,~indx], MPCdmSgrad[x,1,indx]) for x in range(0,MPCdmSgrad[:,0,0].shape[0]) ])
G3_cd = np.array([cohen_d(MPCdmSgrad[x,2,~indx], MPCdmSgrad[x,2,indx]) for x in range(0,MPCdmSgrad[:,0,0].shape[0]) ])

Gcd = [G1_cd, G2_cd, G3_cd]
# plot the gradients
labels=['G'+str(x)+'-cD' for x in list(range(1,Ngrad+1))]
plot_hemispheres(inf_lh, inf_rh, array_name=Gcd, cmap='vlag', nan_color=(0, 0, 0, 1),
  zoom=1.3, size=(900, 750), embed_nb=False, color_range=(-1,1),
  color_bar='right', 
  label_text={'left': labels},
  screenshot=False, filename='/home/bic/rcruces/Desktop/MPC_G3_ASD-HC_cohenD.png')


plot_hemispheres(inf_lh, inf_rh, array_name=G1_cd, 
                 size=(900, 750), zoom=1.3, 
                 embed_nb=False, interactive=False, share='both',
                 nan_color=(0, 0, 0, 1), cmap='vlag', transparent_bg=True, 
                 label_text={'left':['G1-cohenD']}, color_range=(-1, 1),
                 layout_style='grid', color_bar='right',
                 screenshot=True, filename='/home/bic/rcruces/Desktop/cihr_cohend.png', scale=3)



# -----------------------------------------------------------------------------
# ALL subjects G1 BlRd_u
G1 = np.mean(MPCdm[:,0,:], axis=1)
plot_hemispheres(inf_lh, inf_rh, array_name=map_to_labels5k(G1, mask_surf), cmap='RdBu_r',
                  embed_nb=False,  label_text={'left':['G1']}, color_bar='right',
                  layout_style='grid', color_range=(-0.13, 0.13),
                  zoom=1.3, size=(900, 750), nan_color=(0, 0, 0, 1), 
  screenshot=True, filename='/home/bic/rcruces/Desktop/MPC_all_G1.png')


# -----------------------------------------------------------------------------
# t-Test
# Reorder grandientsin three matrices G={subject x Gn}
Nsubj = MPCdm[:,0,:].shape[1]
G1s = np.empty((Nsubj, 9684), float)
for i in range(Nsubj):
    G1s[i,:] = map_to_labels5k(MPCdm[:,0,i], mask_surf)



data_tb = {'sub':[x.split('_')[0].split('sub-')[1] for x in ids],
     'group':[x.split('_')[0].split('sub-')[1][2] for x in ids],
     'ses':[x.split('_')[1].split('ses-')[1] for x in ids],
    'id':ids }
data_df = pd.DataFrame(data_tb)

# terms
term_grp = FixedEffect(data_df['group'])

# contrast
# 1: control, 2: patient
contrast_grp = (data_df.group == 'A').astype(int) - (data_df.group == 'C').astype(int)

# Model
mod = term_grp
Surf = combine_surfaces(inf_lh, inf_rh)

# fitting the model
slm = SLM(
    mod,
    contrast_grp,
    mask=mask_surf,
    surf=Surf,
    correction=['rft'],
    two_tailed=True,
    cluster_threshold=0.05
)
slm.fit(G1s)

# Plot T-values
Save=True
plot_hemispheres(inf_lh, inf_rh, array_name=slm.t[0]*mask_surf, 
                 size=(900, 750), zoom=1.3, 
                 embed_nb=False, interactive=False, share='both',
                 nan_color=(0, 0, 0, 1), cmap='bwr', transparent_bg=True, 
                 label_text={'left':['G1-tval']}, color_range=(-2, 2),
                 layout_style='grid', color_bar='right',
                 screenshot=Save, filename='/home/bic/rcruces/Desktop/cihr_tvals.png', scale=3)

# Plot cluster p-values
Thr = 0.05
plot_pvalues = [np.copy(slm.P["pval"]["C"])]
[np.place(x, np.logical_or(x > Thr, ~mask_surf), np.nan) for x in plot_pvalues]
plot_hemispheres(inf_lh, inf_rh, array_name=plot_pvalues,
                size=(900, 750), zoom=1.3, 
                embed_nb=False, interactive=False, share='both',
                nan_color=(0, 0, 0, 1), cmap='plasma_r', transparent_bg=True, 
                label_text={'left':['G1-pval']}, color_range=(0, Thr),
                layout_style='grid', color_bar='right',
                screenshot=Save, filename='/home/bic/rcruces/Desktop/cihr_pvals.png', scale=3)


# -----------------------------------------------------------------------------
# Group Difference G1 (Cohen's D or T-val?)

#------------------------------------------------------------------------------
# Plot each subject's G1
for i in range(MPCdm[:,0,:].shape[1]):
    plot_hemispheres(inf_lh, inf_rh, array_name=map_to_labels5k(MPCdm[:,0,i], mask_surf), cmap='Spectral_r',
                      embed_nb=False,  label_text={'left':[ids[i]]}, color_bar='right',
                      layout_style='grid', color_range=(-0.13, 0.13),
                      zoom=1.3, size=(900, 750), nan_color=(0, 0, 0, 1), 
      screenshot=False, filename='/home/bic/rcruces/Desktop/MPC_G1_'+ids[i]+'.png')