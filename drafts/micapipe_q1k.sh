#!/bin/bash
#
# Script to run the micapipe Q1K on BIC-SGE
# for i in `ls sub*/ses*/QC/*proc_flair*`; do
# sub=$(echo $i | awk -F '/' '{print $1}')
# ses=$(echo $i | awk -F '/' '{print $2}')
# logs=/data_/mica2/tmpDir/${sub}_${ses}
# qsub -q mica.q -pe smp 10 -l h_vmem=6G -N ${sub}${ses} -e ${logs}.e -o ${logs}.txt /host/yeatman/local_raid/rcruces/git_here/preproc_reports/drafts/micapipe_MICs.sh $sub $ses
# done

sub=$1
ses=$2

# micapipe version
version=v0.2.3

# 2.  Singularity image
img_singularity=/data/mica1/01_programs/micapipe-v0.2.0/micapipe_"${version}".sif

# 3. micapipe command
# Local variables
bids=/data_/mica3/BIDS_Q1K/rawdata
fs_lic=/data_/mica1/01_programs/freesurfer-7.3.2/license.txt
out=/data_/mica3/BIDS_MICs/derivatives
threads=30
tmpDir=/data/mica2/temporaryNetworkProcessing

# Create command string
command="singularity run --writable-tmpfs --containall -B ${bids}:/bids -B ${out}:/out -B ${tmpDir}:/tmp -B ${fs_lic}:/opt/licence.txt ${img_singularity}"

# Run pipeline
${command} \
-bids /bids -out /out -fs_licence /opt/licence.txt -threads ${threads} -sub ${sub} -ses ${ses} \
-proc_structural -proc_surf -post_structural \
-proc_dwi \
-dwi_main /bids/${sub}/${ses}/dwi/${sub}_${ses}_acq-multib38_dir-AP_dwi.nii.gz,/bids/${sub}/${ses}/dwi/${sub}_${ses}_acq-multib70_dir-AP_dwi.nii.gz -regSynth \
-GD -proc_func \
-mainScanStr task-cross_echo-1_bold,task-cross_echo-2_bold,task-cross_echo-3_bold \
-func_pe /bids/${sub}/${ses}/fmap/${sub}_${ses}_dir-AP_epi.nii.gz \
-func_rpe /bids/${sub}/${ses}/fmap/${sub}_${ses}_dir-AP_epi.nii.gz \
-SC -tracts 40M

# Run task-cloudy fmri
${command} \
-bids /bids -out /out -fs_licence /opt/licence.txt -threads ${threads} -sub ${sub} -ses ${ses} \
-proc_func \
-mainScanStr task-cloudy_echo-1_bold,task-cloudy_echo-2_bold,task-cloudy_echo-3_bold \
-func_pe /bids/${sub}/${ses}/fmap/${sub}_${ses}_dir-AP_epi.nii.gz \
-func_rpe /bids/${sub}/${ses}/fmap/${sub}_${ses}_dir-AP_epi.nii.gz \
-QC_subj -QC