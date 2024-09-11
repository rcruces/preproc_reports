#!/bin/bash
#
# HCP - SC with singularity
#

# Subject
sub=$1

# Threads
threads=$2

# BIDS directory
bids=/data/mica3/BIDS_HCP/rawdata

# Out directory
out=/data/mica3/BIDS_HCP/derivatives

# Freesurfer license
fs_lic=/data_/mica1/01_programs/freesurfer-7.3.2/license.txt

# Singularity Image
sif=${MICAPIPE}/micapipe_v0.2.3.sif

# Temporary directory
tmpDir=/host/bb-compx-03/export02/tmp

# Run singularity
singularity run --writable-tmpfs --cleanenv -B ${bids}:/bids -B ${out}:/out -B ${tmpDir}:/tmp -B ${fs_lic}:/opt/license.txt ${sif} -bids /bids -out /out -fs_licence /opt/license.txt -threads 150 -sub ${sub} -SC
