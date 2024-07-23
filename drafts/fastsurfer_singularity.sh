#!/bin/bash
#
# Fastsurfer singularity Example
#

# Image path
fastsurfer_img=/data_/mica1/01_programs/fastsurfer/fastsurfer-cpu-v2.2.0.sif

# Temporary directory path
export TMPDIR=/host/yeatman/local_raid/rcruces/tmp/tmpfiles
tmp=${TMPDIR}

# Path to outputs
SUBJECTS_DIR=/host/yeatman/local_raid/rcruces/tmp/fastsurfer

# Subject ID
idBIDS=MNI152_T1_1mm

# Freesurfer license
fs_license=/data_/mica1/01_programs/freesurfer-7.3.2/license.txt

# MRI|IMG to process
t1=/host/yeatman/local_raid/rcruces/tmp/MNI152_T1_1mm.nii.gz

# Number of threads
threads=15

singularity exec --writable-tmpfs --containall \
                      -B "${SUBJECTS_DIR}":/output \
                      -B "${tmp}":/tmpdir \
                      -B "${t1}":/tmpdir/${idBIDS}_T1w.nii.gz \
                      -B "${fs_license}":/output/license.txt \
                      "${fastsurfer_img}" \
                      /fastsurfer/run_fastsurfer.sh \
                      --fs_license /output/license.txt \
                      --t1 /tmpdir/${idBIDS}_T1w.nii.gz \
                      --sid "${idBIDS}" --sd /output --no_fs_T1 \
                      --parallel --threads "${threads}"
