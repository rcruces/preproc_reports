#!/bin/bash
# ------------------------------------------------
# CHUJS processing workflow
# ------------------------------------------------

# ------------------------------------------------
# Convert DICOM to NIFTI

# ------------------------------------------------
# Set the command to run micapipe
sub=001
ses=01

# micapipe version
version=v0.2.3

# Path to singularity image
img_singularity=/programs/micapipe-v0.2.0/micapipe_"${version}".sif

# BIDS dataset
bids=/P42_T1w/bids_dataset

# freesurfer license
fs_lic=/programs/freesurfer-7.3.2/license.txt

# Output directory (derivatives)
mkdir /P42_T1w/derivatives
out=/P42_T1w/derivatives

# Temporary directory
tmpDir=/tmpDir

# Number of threads
threads=10

# Command to run micapipe
command="singularity run --writable-tmpfs --containall -B ${bids}:/bids -B ${out}:/out -B ${tmpDir}:/tmp -B ${fs_lic}:/opt/licence.txt ${img_singularity}"

# Run micapipe v0.2.3 with singularity
# Structural processing, 
# surface processing, 
# post structural processing,
# Geodesic distance and microstructural processing
# Subject Quality Control reports
${command} \
-bids /bids -out /out -fs_licence /opt/licence.txt -threads ${threads} \
-sub ${sub} -ses ${ses} \
-proc_structural -uni -T1wStr UNIT1,inv-1_MP2RAGE,inv-2_MP2RAGE -mf 12 \
-proc_surf -post_structural -GD \
-MPC -microstructural_img /bids/sub-${sub}/ses-${ses}/anat/sub-${sub}_ses-${ses}_T1map.nii  \
-microstructural_reg FALSE -mpc_acq T1map \
-QC_subj

