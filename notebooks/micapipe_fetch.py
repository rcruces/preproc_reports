import os
import glob
import requests
import tempfile
import nibabel as nib
import numpy as np
import pandas as pd
import seaborn as sns
from brainspace.plotting.surface_plotting import plot_surf
from typing import List, Union, Tuple

from brainspace.mesh.mesh_io import read_surface
from brainspace.plotting import plot_hemispheres

from brainstat.datasets import fetch_mask, fetch_template_surface

import matplotlib as mpl
import matplotlib.pyplot as plt

def fetch_surface(surf_name='fsLR-5k.L.inflated.surf.gii', is_surf=True, nibabel=False):

        """
        Fetches a GIFTI surface file from the micapipe GitHub repository.

        Args:
            surf_name (str): The name of the surface file to download (default is 'fsLR-5k.L.inflated.surf.gii').

        Returns:
            nibabel.Giftidata: The loaded surface data from the GIFTI file.
        """

        # Construct the URL for the surface file on GitHub
        url = f'https://raw.githubusercontent.com/MICA-MNI/micapipe/refs/heads/master/surfaces/{surf_name}'

        # Step 1: Download the surface file from the URL
        response = requests.get(url)

        # Ensure the request was successful
        if response.status_code != 200:
            raise Exception(f"Failed to download the surface file: {surf_name} (Status code: {response.status_code})")

        # Step 2: Save the downloaded file in a temporary directory
        with tempfile.NamedTemporaryFile(delete=False, suffix='.surf.gii') as temp_file:
            temp_file.write(response.content)
            temp_file_name = temp_file.name  # Get the temporary file name

        # Step 3: Read the surface data from the downloaded file (assuming GIFTI format)
        if is_surf:
            if nibabel:
                surf_data = nib.load(temp_file_name)
            else:
                surf_data = read_surface(temp_file_name, itype='gii')
        else:
            surf_data = nib.load(temp_file_name).darrays[0].data

        # Step 4: Remove the temporary file after reading (to avoid cluttering disk)
        os.remove(temp_file_name)

        # Return the surface data
        return surf_data
    
def load_qmri_group(qmri='', subjects_list=None, maps='maps'):
    '''
    This function loads and plots the qMRI intensity maps on fsLR32k midthickness
    '''
    # List the files
    files_lh = []
    for path in subjects_list:
        files_lh.extend(glob.glob(f'{path}/{maps}/*_hemi-L_surf-fsLR-32k{qmri}.func.gii'))
    files_rh = []
    for path in subjects_list:
        files_rh.extend(glob.glob(f'{path}/{maps}/*_hemi-R_surf-fsLR-32k{qmri}.func.gii'))    
    
    files_lh = sorted(files_lh)
    files_rh = sorted(files_rh)

    # Load all the thickness data
    Nqmri=np.concatenate((nib.load(files_lh[0]).darrays[0].data, nib.load(files_rh[0]).darrays[0].data), axis=0).shape[0]

    surf_map=np.empty([len(files_lh), Nqmri], dtype=float)
    for i, _ in enumerate(files_lh):
        surf_map[i,:] = np.hstack(np.concatenate((nib.load(files_lh[i]).darrays[0].data, nib.load(files_rh[i]).darrays[0].data), axis=0))

    # Mean matrix across the x axis (vertices)
    map_mean = np.mean(surf_map, axis=0)
    
    N = surf_map.shape[0]
    print(f"Numer of {qmri} maps: {N}")               
    return(map_mean,surf_map)

class micapipeName:
    ALLOWED_DIRS = ["anat", "dist", "dwi", "func", "maps", "surf", "parc"]

    def __init__(self, **kwargs):
        self.keys = [
            "sub", "ses", "hemi", "space", "surf", "from", "to", "label", 
            "smooth", "desc", "feat", "Dir"
        ]
        self.values = kwargs

    def build(self):
        parts = []
        dir_part = ""
        feat_part = ""
        for key in self.keys:
            if key in self.values:
                if key == "Dir":
                    suffix = self.values[key]
                    if suffix != "Dir":
                        dir_part = suffix
                elif key == "feat":
                    feat_part = f"_{self.values[key]}"
                elif isinstance(self.values[key], (str, int)):
                    parts.append(f"{key}-{self.values[key]}")
        return f"sub-{self.values['sub']}/ses-{self.values['ses']}/{dir_part}/" + "_".join(parts) + feat_part


class micapipe_fetch:
    def __init__(self, derivatives, sub, ses):
        self.derivatives = derivatives
        self.sub = sub
        self.ses = ses

    def surf(self, hemi="both", space="nativepro", surf="fsnative", label="pial", nibabel=False):
        if hemi not in ["L", "R", "both"]:
            raise ValueError("Invalid hemi value. Choose from 'L', 'R', or 'both'.")

        surfaces = {}

        if hemi in ["L", "both"]:
            file_str_lh = micapipeName(sub=self.sub, ses=self.ses, Dir="surf", hemi="L", space=space, surf=surf, label=label).build()
            
            if nibabel==False:
                surfaces["L"] = read_surface(f"{self.derivatives}/{file_str_lh}.surf.gii", itype='gii')
            else:
                surfaces["L"] = nib.load(f"{self.derivatives}/{file_str_lh}.surf.gii").darrays[0].data
                
        if hemi in ["R", "both"]:
            file_str_rh = micapipeName(sub=self.sub, ses=self.ses, Dir="surf", hemi="R", space=space, surf=surf, label=label).build()
            if nibabel==False:
                surfaces["R"] = read_surface(f"{self.derivatives}/{file_str_rh}.surf.gii", itype='gii')
            else:
                surfaces["R"] = nib.load(f"{self.derivatives}/{file_str_rh}.surf.gii").darrays[0].data

        return surfaces

    def feat(self, hemi="both", surf="fsnative", label="pial", feat=None, concatenate=False):
        if hemi not in ["L", "R", "both"]:
            raise ValueError("Invalid hemi value. Choose from 'L', 'R', or 'both'.")

        if feat is None:
            raise ValueError("Feature (feat) must be specified.")

        features = {}

        for hemisphere in ["L", "R"]:
            if hemi in [hemisphere, "both"]:
                kwargs = {
                    "sub": self.sub,
                    "ses": self.ses,
                    "Dir": "maps",
                    "hemi": hemisphere,
                    "surf": surf,
                    "label": "thickness" if feat == "thickness" else label
                }
                if feat != "thickness":
                    kwargs["feat"] = feat

                file_str = micapipeName(**kwargs).build()
                features[hemisphere] = nib.load(f"{self.derivatives}/{file_str}.func.gii").darrays[0].data

        if concatenate:
            feat_cat = np.concatenate((features.get("L", []), features.get("R", [])), axis=0)
            return feat_cat
        else:
            return features


# Functions for plotting
def plot_surfs(
    surfaces,
    values: List[np.ndarray],
    views: Union[List[str], None] = None,
    size: Union[int, Tuple[int, int], None] = None,
    zoom: Union[float, List[float]] = 1.75,
    color_bar="bottom",
    share="both",
    color_range=(-2, 2),
    cmap="rocket",
    transparent_bg=False,
    **kwargs,
):
    """
    surfaces = [hip_mid_l, hip_unf_l, hip_unf_r, hip_mid_r]  Can be 1 or more
    views = ['dorsal', 'lateral', 'lateral', 'dorsal'] Can be 1 or more
    """

    # Append values to surfaces
    my_surfs = {}
    array_names = []
    for i, surf in enumerate(surfaces):
        name = f"surf{i + 1}"
        surf.append_array(values[i], name=name)

        my_surfs[name] = surf
        array_names.append(name)

    # Set the layout as the list of keys from my_surfs
    layout = [list(my_surfs.keys())]

    if size is None:
        size = (200 * len(surfaces), 350)

    return plot_surf(
        my_surfs,
        layout,
        array_name=array_names,
        view=views,
        color_bar=color_bar,
        color_range=color_range,
        share=share,
        cmap=cmap,
        zoom=zoom,
        size=size,
        transparent_bg=transparent_bg,
        **kwargs,
    )

def ridgeplot(matrix, matrix_df=None, Cmap='Reds', color_range=None, 
              Xlab="values", save_path=None, title=None, mean_line=False,
              alpha=0.25, size=(2.5, 3.5), label_pos=0.05, label_col='black'):
    """
    Plot a ridgeplot of the given matrix.

    Parameters:
        matrix (numpy array): The matrix to be plotted.
        matrix_df (pandas DataFrame, optional): DataFrame with additional information about the labels 'id'.
        Cmap (str, optional): Colormap for the ridgeplot. Default is 'Reds'.
        color_range (tuple, optional): Range of values for the x-axis. Default is (matrix.min(), matrix.max()).
        Xlab (str, optional): Label for the x-axis. Default is "values".
        save_path (str, optional): File path to save the plot. If None, the plot is displayed.
        title (str, optional): Title of the plot.
        mean_line (bool, optional): Whether to plot the mean line for all distributions. Default is False.
        alpha (float, optional): Transparency for the KDE plot. Default is 0.25.
        size (tuple, optional): Size of each subplot. Default is (2.5, 3.5).

    Returns:
        None
    """
    # Validate inputs
    if not isinstance(matrix, np.ndarray):
        raise ValueError("matrix must be a numpy array.")
    if matrix_df is not None and not isinstance(matrix_df, pd.DataFrame):
        raise ValueError("matrix_df must be a pandas DataFrame or None.")

    # Create a default DataFrame if none is provided
    if matrix_df is None:
        matrix_df = pd.DataFrame({'id': [f'{i+1}' for i in range(matrix.shape[0])]})

    # Set color range
    if color_range is None:
        color_range = (matrix.min(), matrix.max())

    # Sort matrix rows by their mean values
    mean_row_values = np.mean(matrix, axis=1)
    sorted_indices = np.argsort(mean_row_values)
    sorted_matrix = matrix[sorted_indices]
    sorted_id_x = matrix_df['id'].values[sorted_indices]

    # Prepare data for plotting
    subject = np.repeat(np.arange(1, sorted_matrix.shape[0] + 1), sorted_matrix.shape[1])
    id_x = np.repeat(sorted_id_x, sorted_matrix.shape[1])
    ai = sorted_matrix.flatten()

    df = pd.DataFrame({'feature': ai, 'subject': subject, 'id_x': id_x})

    # Create subplots
    n_rows = sorted_matrix.shape[0]
    fig, axs = plt.subplots(nrows=n_rows, figsize=(3.468504 * size[0], 2.220472 * size[1]), sharex=True, sharey=True)
    fig.set_facecolor('none')

    x = np.linspace(color_range[0], color_range[1], 100)

    for i, ax in enumerate(axs, 1):
        # KDE plot
        sns.kdeplot(df[df["subject"] == i]['feature'], fill=True, color="w", alpha=alpha, linewidth=1.5, legend=False, ax=ax)

        # Set x-axis limits
        ax.set_xlim(color_range[0], color_range[1])

        # Add colormap overlay
        im = ax.imshow(np.vstack([x, x]), cmap=Cmap, aspect="auto", extent=[*ax.get_xlim(), *ax.get_ylim()])
        path = ax.collections[0].get_paths()[0]
        patch = mpl.patches.PathPatch(path, transform=ax.transData)
        im.set_clip_path(patch)

        # Remove spines
        ax.spines[['left', 'right', 'bottom', 'top']].set_visible(False)

        # Adjust ticks and labels
        if i != n_rows:
            ax.tick_params(axis="x", length=0)
        else:
            ax.set_xlabel(Xlab)
        ax.set_yticks([])
        ax.set_ylabel("")
        ax.axhline(0, color="black")
        ax.set_facecolor("none")

        # Add subject label
        ax.text(label_pos, 0.01, sorted_id_x[i - 1], transform=ax.transAxes, fontsize=10, color=label_col, ha='left', va='bottom')

    # Add mean line if requested
    if mean_line:
        mean_asym_all = np.mean(sorted_matrix)
        for ax in axs:
            ax.axvline(x=mean_asym_all, linestyle='dashed', color='black', label=f"Mean: {mean_asym_all:.2f}")

    # Adjust layout
    plt.subplots_adjust(hspace=-0.8)

    # Add title
    if title:
        plt.suptitle(title, y=0.99, fontsize=16)

    # Save or show the plot
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=300)
    else:
        plt.show()