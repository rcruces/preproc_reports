#!/bin/bash
#
# Script to run the micapipe EpiC on BIC-SGE
# for i in `ls sub*/ses*`; do
# sub=$(echo $i | awk -F '/' '{print $1}')
# ses=$(echo $i | awk -F '/' '{print $2}')
# logs=/data_/mica2/tmpDir/${sub}_${ses}
# qsub -q mica.q -pe smp 10 -l h_vmem=6G -N ${sub}${ses} -e ${logs}.e -o ${logs}.txt /host/yeatman/local_raid/rcruces/git_here/preproc_reports/drafts/micapipe_EpiC.sh $sub $ses
# done

sub=$1
ses=$2

# micapipe version
version=v0.2.3

# 2.  Singularity image
img_singularity=/data/mica1/01_programs/micapipe-v0.2.0/micapipe_"${version}".sif

# 3. micapipe command
# Local variables
bids=/data_/mica3/BIDS_EpiC/rawdata
fs_lic=/data_/mica1/01_programs/freesurfer-7.3.2/license.txt
out=/data_/mica3/BIDS_EpiC/derivatives
threads=10
tmpDir=/data_/mica2/tmpDir

# Create command string
command="singularity run --writable-tmpfs --containall -B ${bids}:/bids -B ${out}:/out -B ${tmpDir}:/tmpdir -B ${fs_lic}:/opt/licence.txt ${img_singularity}"

# ------------------------------------------------------------------
# proc_dwi & SC

# cleanup
${command} \
-bids /bids -out /out -fs_licence /opt/licence.txt -threads ${threads} -sub ${sub} -ses ${ses} \
-proc_dwi -SC -acqStr acq-b2000 -cleanup

# reprocess
${command} \
-bids /bids -out /out -fs_licence /opt/licence.txt -threads ${threads} -sub ${sub} -ses ${ses} \
-proc_dwi -SC -dwi_acq b2000 -tmpDir /tmpdir -regSynth -dwi_upsample \
-dwi_main /bids/sub-${sub}/ses-${ses}/dwi/sub-${sub}_ses-${ses}_acq-b2000_dir-AP_dwi.nii.gz


# ------------------------------------------------------------------
# All modules
# ${command} \
# -bids /bids -out /out -fs_licence /opt/licence.txt -threads ${threads} -sub ${sub} -ses ${ses} \
    # -proc_surf -surf_dir /out/freesurfer/${sub}_${ses} \
    # -post_structural -freesurfer \
    # -proc_flair \
    # -GD \
    # -regSynth \
    # -proc_dwi \
    # -SC

# ------------------------------------------------------------------
# func cleanup
# ${command} \
# -bids /bids -out /out -fs_licence /opt/licence.txt -threads ${threads} -sub ${sub} -ses ${ses} \
# -proc_func -acqStr desc-se_task-rest_bold -cleanup

# ${command} \
# -bids /bids -out /out -fs_licence /opt/licence.txt -threads ${threads} -sub ${sub} -ses ${ses} \
# -proc_func -acqStr desc-se_task-sternberg_bold -cleanup

# Run func
# ${command} \
# -bids /bids -out /out -fs_licence /opt/licence.txt -threads ${threads} -sub ${sub} -ses ${ses} \
# -proc_func -mainScanStr task-rest_bold -tmpDir ${tmpDir} -noFIX -NSR

# ${command} \
# -bids /bids -out /out -fs_licence /opt/licence.txt -threads ${threads} -sub ${sub} -ses ${ses} \
# -proc_func -mainScanStr task-sternberg_bold -tmpDir ${tmpDir} -noFIX -NSR
