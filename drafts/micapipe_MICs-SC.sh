#!/bin/bash
#
# Script to run the micapipe MICs on BIC-SGE

# /host/verges/tank/data/tmp | subjects=(PX173 PX152 PX198 PX216 PX050 PX043 PX064 PX123 PX122 PX132)
# /host/bb-compx-03/export02/tmp | subjects=(HC082 HC081 HC135 HC083 HC118 HC076 HC062 HC088 HC130 HC152)
# /data_/mica3/tmp_proc | subjects=(PNC003 PNC006 PNC009 PNC011 PNC016 PNC018 PNC019 PNC024 PNC026 PNC039 PNE006 PNE014 PNE018 PNE021 PNE029 PNE031 PNE041 PNE049 PNE054 PNE055)

# for sub in ${subjects[@]}; do
#   for ses in 01 02; do
#   logs=/host/verges/tank/data/tmp/2026_MICs-SC_${sub}_${ses}
#   qsub -q mica.q -pe smp 10 -l h_vmem=6G -N ${sub}${ses} -e ${logs}.e -o ${logs}.txt /host/yeatman/local_raid/rcruces/git_here/preproc_reports/drafts/micapipe_MICs-SC.sh $sub
#   done
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
threads=10
tmpDir=/data_/mica3/tmp_proc

# Create command string
command="singularity run --writable-tmpfs --containall -B ${bids}:/bids -B ${out}:/out -B ${tmpDir}:/tmpdir -B ${fs_lic}:/opt/licence.txt ${img_singularity}"

echo "--------------------------------------------------------"
echo "tmpDir:   ${tmpDir}"
echo "${command} -bids /bids -out /out -fs_licence /opt/licence.txt -threads ${threads} -sub ${sub} -ses ${ses} -tmpDir /tmpdir -SC"
echo "--------------------------------------------------------"

# SC processing
${command} -bids /bids -out /out -fs_licence /opt/licence.txt -threads ${threads} -sub ${sub} -ses ${ses} -tmpDir /tmpdir -SC -tracts 5M

# SC processing
${command} -bids /bids -out /out -fs_licence /opt/licence.txt -threads ${threads} -sub ${sub} -ses ${ses} -tmpDir /tmpdir -SC -tracts 10M

# SC processing
${command} -bids /bids -out /out -fs_licence /opt/licence.txt -threads ${threads} -sub ${sub} -ses ${ses} -tmpDir /tmpdir -SC -tracts 20M


