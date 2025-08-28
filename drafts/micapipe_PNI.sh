#!/bin/bash
#
# Script to run the A subject with singularity on BIC-SGE
# Database: Presicion Neuroimaging 7T
#
# ----------------------------------------------------------

# Subjects and sessions
sub=$1
ses=$2

# micapipe version
version=v0.2.3

# 2.  Singularity image
img_singularity=/data/mica1/01_programs/micapipe-v0.2.0/micapipe_"${version}".sif

# 3. micapipe command
# Local variables
bids=/data_/mica3/BIDS_PNI/rawdata
fs_lic=/data_/mica1/01_programs/freesurfer-7.3.2/license.txt
out=/data_/mica3/BIDS_PNI/derivatives
tmpDir=/data/mica2/tmpDir
threads=15

# Create command string
command="singularity run --writable-tmpfs --containall -B ${bids}:/bids -B ${out}:/out -B ${tmpDir}:/tmp -B ${fs_lic}:/opt/licence.txt ${img_singularity}"

# MPC cleanup
${command} \
-bids /bids -out /out -fs_licence /opt/licence.txt -threads ${threads} -sub ${sub} -ses ${ses} \
-MPC -cleanup -acqStr T2starmap_proc
