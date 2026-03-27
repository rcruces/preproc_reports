#!/bin/bash
#
# Script to run the micapipe deface algorithm on Healthy PET mk6240
# for i in sub*/ses*; do
# sub=$(echo ${i/sub-} | awk -F '/' '{print $1}')
# ses=$(echo ${i/ses-} | awk -F '/' '{print $2}')
# micapipe_anonymize -sub ${sub} -ses ${ses} -out ${out} -bids ${bids} -deface -robust -threads ${threads}
# done

sub=$1
ses=$2

# micapipe version
version=v0.2.3

# 2.  Singularity image
img_singularity=/data/mica1/01_programs/micapipe-v0.2.0/micapipe_"${version}".sif

# 3. micapipe command
# Local variables
bids=/data_/mica3/MICA-PET/BIDS_MK6240_HC/rawdata
fs_lic=/data_/mica1/01_programs/freesurfer-7.3.2/license.txt
out=/data_/mica3/MICA-PET/BIDS_MK6240_HC/derivatives
threads=15
tmpDir=/tmp

# Create command string
command="singularity exec --writable-tmpfs --containall -B ${bids}:/bids -B ${out}:/out -B ${tmpDir}:/tmp -B ${fs_lic}:/opt/licence.txt ${img_singularity}"

# Run Command
${command} /opt/micapipe/functions/micapipe_anonymize -sub ${sub} -ses ${ses} -out ${out} -bids ${bids} -deface -robust -threads ${threads}