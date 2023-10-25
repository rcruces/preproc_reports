#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 15:30:07 2023

@author: rcruces
"""
import os
import glob
import seaborn as sns
import numpy as np
import nibabel as nb
from brainspace.plotting import plot_hemispheres
from brainspace.mesh.mesh_io import read_surface
import matplotlib.pyplot as plt
import pandas as pd
from brainstat.datasets.base import combine_surfaces

def plot_connectome(mtx, Title='matrix plot', xlab='X', ylab='Y', col='rocket', vmin=None, vmax=None,
                   xticklabels='auto', yticklabels='auto',xrot=90, yrot=0, save_path=None):

    '''
    This optional function, only plots a connectome as a heatmap
    Parameters
    ----------
    mtx : np.array
    Returns
    -------
    '''
    f, ax = plt.subplots(figsize=(15,10))
    g = sns.heatmap(mtx, ax=ax, cmap=col, vmin=vmin, vmax=vmax, xticklabels=xticklabels, yticklabels=yticklabels)
    g.set_xlabel(xlab)
    g.set_ylabel(ylab)
    g.set_title(Title)
    # Rotate the x-axis labels
    # rotate tick labels
    ax.set_xticklabels(ax.get_xticklabels(), rotation=xrot, ha='right')
    ax.set_yticklabels(ax.get_yticklabels(), rotation=yrot, ha='right')

    if save_path is not None:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')

# -----------------------------------------------------------------------------
# plot ridge plot
def plot_ridgeplot(matrix, matrix_df=None, Cmap='rocket', Range=(2.5, 4.5), Xlab="flair", save_path=None, title=None):
    import numpy as np
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    mean_row_values = np.mean(matrix, axis=1)
    sorted_indices = np.argsort(mean_row_values)
    sorted_matrix = matrix[sorted_indices]
    sorted_id_x = np.array(matrix_df)[sorted_indices]

    ai = sorted_matrix.flatten()
    subject = np.array([])
    id_x = np.array([])

    for i in range(sorted_matrix.shape[0]):
        label = np.array([str(i+1) for j in range(sorted_matrix.shape[1])])
        subject = np.concatenate((subject, label))
        id_label = np.array([sorted_id_x[i] for j in range(sorted_matrix.shape[1])])
        id_x = np.concatenate((id_x, id_label))

    d = {'feature': ai,
         'subject': subject,
         'id_x': id_x
        }
    df = pd.DataFrame(d)

    f, axs = plt.subplots(nrows=sorted_matrix.shape[0], figsize=(3.468504*2.5, 2.220472*3.5), sharex=True, sharey=True)
    f.set_facecolor('none')

    x = np.linspace(Range[0], Range[1], 100)

    for i, ax in enumerate(axs, 1):
        sns.kdeplot(df[df["subject"]==str(i)]['feature'],
                    fill=True,
                    color="w",
                    alpha=0.25,
                    linewidth=1.5,
                    legend=False,
                    ax=ax)
        
        ax.set_xlim(Range[0], Range[1])
        
        im = ax.imshow(np.vstack([x, x]),
                       cmap=Cmap,
                       aspect="auto",
                       extent=[*ax.get_xlim(), *ax.get_ylim()]
                      )
        ax.collections
        path = ax.collections[0].get_paths()[0]
        patch = mpl.patches.PathPatch(path, transform=ax.transData)
        im.set_clip_path(patch)
           
        ax.spines[['left','right','bottom','top']].set_visible(False)
        
        if i != sorted_matrix.shape[0]:
            ax.tick_params(axis="x", length=0)
        else:
            ax.set_xlabel(Xlab)
            
        ax.set_yticks([])
        ax.set_ylabel("")
        
        ax.axhline(0, color="black")

        ax.set_facecolor("none")

    for i, ax in enumerate(axs):
        if i == sorted_matrix.shape[0] - 1:
            ax.set_xticks([Range[0], Range[1]])  # Set x-axis ticks for the bottom plot
        else:
            ax.set_xticks([])  # Remove x-axis ticks from other plots
        ax.text(0.05, 0.01, sorted_id_x[i], transform=ax.transAxes, fontsize=10, color='black', ha='left', va='bottom')

    plt.subplots_adjust(hspace=-0.8)
    
    if title:
        plt.suptitle(title, y=0.99, fontsize=16)

    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=300)
    else:
        plt.show()
        
def NormalizeData(data):
    return (data - np.min(data)) / (np.max(data) - np.min(data))

from scipy import stats
def NormalizeMode(data):   
    return(data / stats.mode(data, keepdims=True)[0][0])

os.chdir('/data_/mica3/BIDS_MICs/derivatives/micapipe_v0.2.0')




# Path to MICAPIPE
micapipe=os.popen("echo $MICAPIPE").read()[:-1]
i5_lh = read_surface(micapipe + '/surfaces/fsLR-5k.L.inflated.surf.gii', itype='gii')
i5_rh = read_surface(micapipe + '/surfaces/fsLR-5k.R.inflated.surf.gii', itype='gii')
mask_lh =  nb.load(micapipe + '/surfaces/fsLR-5k.R.mask.shape.gii').darrays[0].data
mask_rh =  nb.load(micapipe + '/surfaces/fsLR-5k.L.mask.shape.gii').darrays[0].data
mask_10k = np.concatenate((mask_lh, mask_rh), axis=0)

# Labels and boolean mask
mask_surf = mask_10k != 0

# Load all flair data
lh_str='sub-PX*/ses-01/maps/*hemi-L_surf-fsLR-5k_label-midthickness_flair*'
rh_str='sub-PX*/ses-01/maps/*hemi-R_surf-fsLR-5k_label-midthickness_flair*'
lh_files=sorted(glob.glob(lh_str))
rh_files=sorted(glob.glob(rh_str))

# get the BIDS ids
bids_ids = [elem.split('/')[0] + '_' + elem.split('/')[1] for elem in lh_files]

# Load all the flair data
Nth=np.concatenate((nb.load(lh_files[0]).darrays[0].data, nb.load(rh_files[0]).darrays[0].data), axis=0).shape[0]

surf_map=np.empty([len(lh_files), Nth], dtype=float)
for i, _ in enumerate(lh_files):
    #print(f)
    surf_map[i,:] = np.hstack(np.concatenate((nb.load(lh_files[i]).darrays[0].data, nb.load(rh_files[i]).darrays[0].data), axis=0))

# Normalize data
surf_norm = np.apply_along_axis(NormalizeMode, 1, surf_map)

# Mean matrix across the x axis (vertices)
map_mean = np.mean(surf_norm, axis=0)

# Plot the mean FEATURE 10mm on conte69 surface
quantile=(0.15, 0.999)
cmap='afmhot'
Range=(np.quantile(map_mean, quantile[0]), np.quantile(map_mean, quantile[1]))

plot_hemispheres(i5_lh, i5_rh, array_name=map_mean, cmap=cmap, nan_color=(0, 0, 0, 1),
                   zoom=1.3, size=(900, 750), embed_nb=False,
                   color_bar='right', layout_style='grid', color_range=Range,
                   label_text={'left': ['Lateral', 'Medial'], 'top': ['Left', 'Right']},
                   screenshot=False)

# Plot the mean VARIANCE 10mm on conte69 surface
map_var = np.var(surf_map, axis=0)

plot_hemispheres(i5_lh, i5_rh, array_name=map_var, cmap='Spectral_r', nan_color=(0, 0, 0, 1),
                     zoom=1.3, size=(900, 750), embed_nb=False,
                     color_bar='right', layout_style='grid', color_range=(0, np.nanquantile(map_var, 0.95)),
                     label_text={'left': ['Lateral', 'Medial'], 'top': ['Left', 'Right']},
                     screenshot=False)

## correlation matrix
corr = np.corrcoef(surf_map)
plot_connectome(corr, 'Subject similarity', xlab=None, ylab=None, col='Spectral_r', vmin=0.1, vmax=1,
                           yticklabels=bids_ids, xticklabels=bids_ids)


plot_ridgeplot(surf_map, bids_ids, title="flair distribution in TLE", Cmap='afmhot', Range=(0,6))




plot_ridgeplot(surf_norm, bids_ids, title="flair distribution in TLE", Cmap='afmhot', Range=(0,1))


# Ttest
px_flair = np.copy(surf_norm)
px_ids = bids_ids
hc_flair = np.copy(surf_norm)
hc_ids = np.copy(bids_ids)

ids = px_ids + list(hc_ids)


flair = np.vstack((px_flair, hc_flair))

data_tb = {'sub':[x.split('_')[0].split('sub-')[1] for x in ids],
     'group':[x.split('_')[0].split('sub-')[1][0] for x in ids],
     'ses':[x.split('_')[1].split('ses-')[1] for x in ids],
    'id':ids }
data_df = pd.DataFrame(data_tb)

from brainstat.stats.SLM import SLM
from brainstat.stats.terms import FixedEffect

# terms
term_grp = FixedEffect(data_df['group'])

# contrast
# 1: control, 2: patient
contrast_grp = (data_df.group == 'P').astype(int) - (data_df.group == 'H').astype(int)

# Model
mod = term_grp
Surf = combine_surfaces(i5_lh, i5_rh)

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
slm.fit(flair)

# Plot T-values
Save=False
plot_hemispheres(i5_lh, i5_rh, array_name=slm.t[0]*mask_surf, 
                 size=(900, 750), zoom=1.3, 
                 embed_nb=False, interactive=False, share='both',
                 nan_color=(0, 0, 0, 1), cmap='vlag', transparent_bg=True, 
                 label_text={'left':['G1-tval']}, color_range='sym',
                 layout_style='grid', color_bar='right',
                 screenshot=Save, filename='/home/bic/rcruces/Desktop/cihr_tvals.png', scale=3)

# Plot cluster p-values
Thr = 0.025
plot_pvalues = [np.copy(slm.P["pval"]["C"])]
[np.place(x, np.logical_or(x > Thr, ~mask_surf), np.nan) for x in plot_pvalues]
plot_hemispheres(i5_lh, i5_rh, array_name=plot_pvalues,
                size=(900, 750), zoom=1.3, 
                embed_nb=False, interactive=False, share='both',
                nan_color=(0, 0, 0, 1), cmap='plasma_r', transparent_bg=True, 
                label_text={'left':['G1-pval']}, color_range=(0, Thr),
                layout_style='grid', color_bar='right',
                screenshot=Save, filename='/home/bic/rcruces/Desktop/cihr_pvals.png', scale=3)

# Index of the controls
indx =np.array(["HC" in f for f in ids])

cD = np.array([cohen_d(flair[~indx,x], flair[indx,x]) for x in range(0,flair.shape[1]) ])

plot_hemispheres(inf_lh, inf_rh, array_name=cD, 
                 size=(900, 750), zoom=1.3, 
                 embed_nb=False, interactive=False, share='both',
                 nan_color=(0, 0, 0, 1), cmap='vlag', transparent_bg=True, 
                 label_text={'left':['G1-cohenD']}, color_range=(-0.7, 0.7),
                 layout_style='grid', color_bar='right',
                 screenshot=False, filename='/home/bic/rcruces/Desktop/cihr_cohend.png', scale=3)