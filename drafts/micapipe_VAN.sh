#!/bin/bash
#
# Script to run the micapipe on BIDS_VAN data
# 
sub=$1
ses=$2

# micapipe version
version=v0.2.3

# 2.  Singularity image
img_singularity=/data/mica1/01_programs/micapipe-v0.2.0/micapipe_"${version}".sif

# 3. micapipe command
# Local variables
bids=/data/mica3/BIDS_VAN/rawdata/
out=/data/mica3/BIDS_VAN/derivatives/
fs_lic=/data_/mica1/01_programs/freesurfer-7.3.2/license.txt
tmpDir=/data/mica2/temporaryNetworkProcessing/
threads=15

# Create command string
command="singularity run --writable-tmpfs --containall -B ${bids}:/bids -B ${out}:/out -B ${tmpDir}:/tmp -B ${fs_lic}:/opt/licence.txt ${img_singularity}"

# All modules
${command} \
-bids /bids -out /out -fs_licence /opt/licence.txt -threads ${threads} -tmpDir /tmp -sub ${sub} -ses ${ses} \
    -proc_surf -surf_dir /out/freesurfer/${sub}_${ses} \
    -post_structural -freesurfer \
    -proc_flair \
    -proc_dwi -regSynth -SC \
    -GD \
    -proc_func -mainScanStr task-rest_bold -NSR -noFIX -dropTR -QC_subj

