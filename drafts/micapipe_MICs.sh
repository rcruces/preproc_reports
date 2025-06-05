#!/bin/bash
#
# Script to run the micapipe MICs on BIC-SGE
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
bids=/data_/mica3/BIDS_MICs/rawdata
fs_lic=/data_/mica1/01_programs/freesurfer-7.3.2/license.txt
out=/data_/mica3/BIDS_MICs/derivatives
threads=15
tmpDir=/tmp

# Create command string
command="singularity run --writable-tmpfs --containall -B ${bids}:/bids -B ${out}:/out -B ${tmpDir}:/tmp -B ${fs_lic}:/opt/licence.txt ${img_singularity}"

# All modules
${command} \
-bids /bids -out /out -fs_licence /opt/licence.txt -threads ${threads} -sub ${sub} -ses ${ses} \
    -proc_surf -surf_dir /out/freesurfer/${sub}_${ses} \
    -post_structural -freesurfer \
    -proc_flair \
    -GD \
    -proc_func -phaseReversalRun 1 -dropTR \
    -MPC -mpc_acq T1map -regSynth \
    -microstructural_reg FALSE \
    -microstructural_img /bids/sub-${sub}/ses-${ses}/anat/*_T1map.nii.gz \
    -proc_dwi \
    -SC

# flair cleanup
# ${command} \
# -bids /bids -out /out -fs_licence /opt/licence.txt -threads ${threads} -sub ${sub} -ses ${ses} \
# -proc_flair -cleanup
