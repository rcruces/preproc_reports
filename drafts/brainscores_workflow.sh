#!/bin/bash
#
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# Dicoms/Niftis to BIDS
# One rawdata BIDS directory withh all subjects or
# Individual directories with each process??
# NAMING convention of the sub-<ID>

# ----------------------------------------------------------------
# T1w enhacement
mri_synthsr --i T1w_cor-reg.nii.gz --o T1w_cor-SynthSR.nii.gz --lowfield --cpu --threads 20

# ----------------------------------------------------------------
# MICAPIPE
# Threads
threads=15

# BIDS directory
bids=/data_/mica3/BIDS_brainscores/rawdata

# Out directory
out=/data_/mica3/BIDS_brainscores/derivatives

# Subject
sub=sub-mx250204

# Freesurfer license
fs_lic=/data_/mica1/01_programs/freesurfer-7.3.2/license.txt

# Singularity Image
sif_micapipe=${MICAPIPE}/micapipe_v0.2.3.sif

# Temporary directory
tmpDir=/data_/mica3/tmp/

# Run singularity
singularity run --writable-tmpfs --cleanenv \
    -B ${bids}:/bids -B ${out}:/out -B ${tmpDir}:/tmp \
    -B ${fs_lic}:/opt/license.txt \
    ${sif_micapipe} -bids /bids -out /out -fs_licence /opt/license.txt \
    -threads ${threads} -sub ${sub} -proc_structural -proc_surf -post_structural -proc_flair -QC_subj

# ----------------------------------------------------------------
# HIPUNFOLD
export HIPPUNFOLD_CACHE_DIR=/host/yeatman/local_raid/rcruces/cache/hippunfold
export SNAKEMAKE_OUTPUT_CACHE=/host/yeatman/local_raid/rcruces/cache/snakemake

# Singularity Image
sif_hipunfold=/data/mica1/01_programs/singularity/hippunfold_v1.4.1.sif

cd $SNAKEMAKE_OUTPUT_CACHE
singularity run \
        -B ${out}/micapipe_v0.2.0:/mica \
        -B ${out}/hippunfold_v1.4.1:/hipp \
        ${sif_hipunfold} /mica /hipp \
        participant \
        --participant_label ${sub/sub-} \
        --modality T1w \
        --filter-T1w space=nativepro \
        --output-density 0p5mm 2mm \
        --cores ${threads} 

# ----------------------------------------------------------------
# ZBRAINS
# Set the path to the dataset, or the folder containing the 'derivatives' folder
dataset_ref="/data/mica3/BIDS_MICs"

# Set the directories for micapipe, hippunfold, and zbrains, which will be looked for in the 'derivates' folder
zbrains_dir="${out}/zbrains"
micapipe_dir="${out}/micapipe_v0.2.0"
hippunfold_dir="${out}/hippunfold_v1.3.0"

# Set the paths to the demographic control and patient files
# The demo_controls are only needed for the analysis, and define the control samples to compare against.
# The demo_patients can be provided if desired, which will run all the patients in the list when the "all" keyword is used,
# otherwise the 'all' keyword will run every patient possible, given the micapipe and hippunfold availability, or, for the analysis
# it will run on all patients who have had zbrains proc run before.
demo_controls="/host/oncilla/local_raid/oualid/zbrains_csvs/participants_mics_hc.csv"

# Set the subject IDs and session IDs to 'all', using all patients defined in the PX_participants file.
sub=sub-mx250204
ses="ses-01"

# The code below runs zbrains preserving the old behaviour, with a smooth_ctx of 10, a smooth_hip of 5, and a label_ctx of 'white'
# The new defaults for this are a smooth_ctx of 5, a smooth_hip of 2, and a label_ctx of 'midthickness'
# Much of the new volumetric code is dependent on cortical midthickness, so it is recommended.
zbrains --run "proc analysis"\
        --sub "${sub}" \
        --ses "${ses}" \
        --dataset "/data_/mica3/BIDS_brainscores" \
        --zbrains ${zbrains_dir} \
        --micapipe ${micapipe_dir} \
        --hippunfold ${hippunfold_dir} \
        --dataset_ref ${dataset_ref} \
        --zbrains_ref "wbrains" \
        --demo_ref ${demo_controls} \
        --column_map participant_id=ID session_id=SES \
        --smooth_ctx 10 \
        --smooth_hip 5 \
        --n_jobs 4 \
        --n_jobs_wb 4 \
        --label_ctx "white" \
        --dicoms 0 \
        --verbose 2 \
        --pyinit=/data/mica1/03_projects/ian/anaconda3 \
	--volumetric 0

# ----------------------------------------------------------------
# MICAFLOW

