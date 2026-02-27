#!/bin/bash
#
# Script to run the micapipe MICs on BIC-SGE


#subjects=(sub-HC001 sub-HC002 sub-HC006 sub-HC007 sub-HC011 sub-HC012 sub-HC013 sub-HC015 sub-HC018 sub-HC020 sub-HC023 sub-HC025 sub-HC026 sub-HC028 sub-HC030 sub-HC037 sub-HC038 sub-HC040 sub-HC041 sub-HC049 sub-HC050 sub-HC054 sub-HC061 sub-HC069 sub-HC070 sub-HC072 sub-HC075 sub-HC077 sub-HC084 sub-HC085 sub-HC088 sub-HC095 sub-HC097 sub-HC098 sub-HC101 sub-HC106 sub-HC107 sub-HC112 sub-HC113 sub-HC115 sub-HC119 sub-HC120 sub-HC121 sub-HC123 sub-HC127 sub-HC130 sub-HC131 sub-HC132 sub-HC133 sub-HC147 sub-HC149 sub-HC152 sub-HC155 sub-HC161 sub-HC166 sub-PX002 sub-PX004 sub-PX008 sub-PX009 sub-PX014 sub-PX015 sub-PX017 sub-PX018 sub-PX019 sub-PX021 sub-PX023 sub-PX026 sub-PX027 sub-PX028 sub-PX030 sub-PX031 sub-PX034 sub-PX038 sub-PX040 sub-PX041 sub-PX043 sub-PX045 sub-PX049 sub-PX051 sub-PX054 sub-PX056 sub-PX061 sub-PX064 sub-PX065 sub-PX066 sub-PX069 sub-PX074 sub-PX077 sub-PX078 sub-PX087 sub-PX094 sub-PX099 sub-PX102 sub-PX104 sub-PX106 sub-PX107 sub-PX121 sub-PX127 sub-PX129 sub-PX133 sub-PX134 sub-PX138 sub-PX140 sub-PX142 sub-PX143 sub-PX145 sub-PX148 sub-PX149 sub-PX150 sub-PX152 sub-PX155 sub-PX156 sub-PX162 sub-PX169 sub-PX170 sub-PX171 sub-PX172 sub-PX174 sub-PX178 sub-PX183 sub-PX184 sub-PX189 sub-PX193 sub-PX195 sub-PX197 sub-PX198 sub-PX209 sub-PX216 sub-PX219 sub-PX224 sub-PX251 sub-PX262 sub-PX263 sub-PX265)

# for i in ${subjects[@]}; do
# sub=$(echo $i | awk -F '/' '{print $1}')
# ses=$(echo $i | awk -F '/' '{print $2}')
# logs=/data_/mica2/tmpDir/2026_MICs-SC_${sub}_${ses}
# qsub -q mica.q -pe smp 10 -l h_vmem=6G -N ${sub}${ses} -e ${logs}.e -o ${logs}.txt /host/yeatman/local_raid/rcruces/git_here/preproc_reports/drafts/micapipe_MICs.sh $sub $ses
# done

sub=$1
ses=01

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

# Create command string
command="singularity run --writable-tmpfs --containall -B ${bids}:/bids -B ${out}:/out -B /export02:/export02 -B /export03:/export03 -B ${TMP}:/tmp -B ${fs_lic}:/opt/licence.txt ${img_singularity}"

# SC
${command} \
-bids /bids -out /out -fs_licence /opt/licence.txt -threads ${threads} -sub ${sub} -ses ${ses} -SC