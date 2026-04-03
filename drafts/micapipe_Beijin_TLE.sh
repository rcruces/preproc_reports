#!/bin/bash
#
#---------------------------------------------------
# For loop to run the micapipe Beijing on BIC-SGE
# for i in sub*/ses*; do
#   sub=$(echo $i | cut -d/ -f1);
#   ses=$(echo $i | cut -d/ -f2);
#   logs=/host/verges/tank/data/tmp/${sub}_${ses}
#   qsub -q mica.q -pe smp 10 -l h_vmem=6G -N ${sub}${ses} -e ${logs}.e -o ${logs}.txt /host/yeatman/local_raid/rcruces/git_here/preproc_reports/drafts/micapipe_Beijing_FCD.sh $sub $ses
# done
#---------------------------------------------------

sub=$1
ses=$2

# micapipe version
version=v0.2.3

# 2.  Singularity image
img_singularity=/data/mica1/01_programs/micapipe-v0.2.0/micapipe_"${version}".sif

# 3. micapipe command
# Dataset paths
bids=/host/verges/tank/data/mo/BIDS_20251201_TLE_BTTH1/rawdata
out=/host/verges/tank/data/mo/BIDS_20251201_TLE_BTTH1/derivatives

# Local variables
fs_lic=/data_/mica1/01_programs/freesurfer-7.3.2/license.txt
threads=10
tmpDir=/host/bb-compx-03/export02/tmp
export TMPDIR=$tmpDir

# Create command string
command="singularity run --writable-tmpfs --containall -B ${bids}:/bids -B ${out}:/out -B ${tmpDir}:/tmpdir -B ${tmpDir}:/tmp -B ${fs_lic}:/opt/licence.txt ${img_singularity}"

# NOt sure if this is needed: echo $tmpDir > ~/.mrtrix.conf 

# Permission (makes sure that the files are rx for aug)
umask 002

echo "--------------------------------------------------------"
echo "tmpDir:   ${tmpDir}"
echo "${command} -bids /bids -out /out -fs_licence /opt/licence.txt -threads ${threads} -sub ${sub} -ses ${ses} -tmpDir /tmpdir"
echo "--------------------------------------------------------"

# ------------------------------------------------------------------
# proc_dwi
${command} \
-bids /bids -out /out -fs_licence /opt/licence.txt -threads ${threads} -sub ${sub} -ses ${ses} \
-proc_dwi -tmpDir /tmpdir -regSynth -dwi_upsample \
-dwi_main /bids/sub-${sub}/ses-${ses}/dwi/sub-${sub}_ses-${ses}_dwi.nii.gz -dwi_rpe ""

# CLEANUP -proc_dwi
# ${command} -bids /bids -out /out -fs_licence /opt/licence.txt -threads ${threads} -sub ${sub} -ses ${ses} \
# -cleanup -proc_dwi