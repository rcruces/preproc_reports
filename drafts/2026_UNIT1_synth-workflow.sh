#!/usr/bin/env bash
# Author: Raul R. Cruces
# ==============================================================================
# Synthetic MRI generation and FastSurfer surface reconstruction workflow
# ==============================================================================
#
# Generate a synthetic T1w-like image from an ultra-high-resolution MP2RAGE
# UNIT1 acquisition and use it for cortical surface reconstruction with
# FastSurfer.
#
# PROCESSING STEPS
# ----------------
# 1. MP2RAGE background denoising
# 2. N4 bias field correction
# 3. Non-local means denoising (parameter sweep)
# 4. SynthSeg segmentation
# 5. SynthSR generation and warping back to native geometry
# 6. FastSurfer cortical reconstruction
# 7. micapipe structural processing
# 8. Intracortical and SWM mapping
#
# EMPIRICAL RESULT
# ----------------
# Tested NLM parameters:
#   sigma = [0,1,2]
#   beta  = [0,1,2]
#
# Visual QC suggested:
#   sigma = 2
#   beta  = 1
#
# provides the best compromise between denoising and preservation of cortical
# detail.
#
# Dependencies:
#   - ANTs
#   - FastSurfer v2.5.4
#   - FreeSurfer
#   - SynthSeg
#   - SynthSR
#   - FSL
#   - Singularity/Apptainer
#   - custom denoiseNLM wrapper


# ------------------------------------------------------------------------------
# PREPROCESSING (performed previously)
# ------------------------------------------------------------------------------
#
# Starting image:  UNIT1
#
# mp2rage_denoising (from micapipe): rec-denoised_acq-05mm_UNIT1.nii.gz
# > Salt-and-pepper background cleanup
#
# N4 bias-field correcrtion: rec-denoisedN4_acq-05mm_UNIT1.nii.gz


# ------------------------------------------------------------------------------
# 3. NON-LOCAL MEANS DENOISING PARAMETER SWEEP
# ------------------------------------------------------------------------------
#
# Evaluate combinations of:
#
#   sigma = noise estimate
#   beta  = smoothing strength
#
# Output naming convention:
#
#   rec-NLMs{sigma}b{beta}
#
# Special case:
#   sigma=0 beta=1 corresponds to the denoiseNLM default settings.
#
# ------------------------------------------------------------------------------
# NOTE: ANTs DenoiseImage was also evaluated but is considerably slower than the

idBIDS=sub-PNE032_ses-a1
export denoiseNLM=/host/yeatman/local_raid/rcruces/git_here/MICA-7t/functions/denoiseNLM

for beta in 0 1 2; do
    for sigma in 0 1 2; do

        if [[ "$sigma" -eq 0 && "$beta" -eq 1 ]]; then
            nlm_entity="rec-NLMs${sigma}b${beta}DEFAULT"
        else
            nlm_entity="rec-NLMs${sigma}b${beta}"
        fi

        denoiseNLM ${idBIDS}_rec-denoisedN4_acq-05mm_UNIT1.nii.gz ${idBIDS}_${nlm_entity}_acq-05mm_UNIT1 . -threads 20 -sigma ${sigma} -beta ${beta}
    done
done

# ------------------------------------------------------------------------------
# SELECT OPTIMAL IMAGE FOR DOWNSTREAM PROCESSING
# ------------------------------------------------------------------------------
# NLM sigma=2 beta=1
# This image is used as the reference space for all subsequent processing.
reg_img=${idBIDS}_rec-NLMs2b1_acq-05mm_UNIT1.nii.gz

# SynthSeg segmentation
reg_img_seg=${idBIDS}_rec-NLMs2b1_acq-1mm_desc-synthSEG_UNIT1.nii.gz

# Resampled volume produced by SynthSeg
reg_img_seg_res=synthseg/${idBIDS}_rec-NLMs2b1_acq-05mm_UNIT1_resampled.nii.gz

# SynthSR output
reg_img_sr=${idBIDS}_rec-NLMs2b1_acq-1mm_desc-synthSR_UNIT1.nii.gz


# ------------------------------------------------------------------------------
# 4. SYNTHSEG: Robust tissue segmentation.
# ------------------------------------------------------------------------------
# SynthSeg additionally produces a resampled image that is later used as the input to SynthSR.

mri_synthseg --i ${reg_img} --o ${reg_img_seg} --robust --threads 20 --cpu --resample synthseg


# ------------------------------------------------------------------------------
# 5. SYNTHSR: SynthSR generates a synthetic (NON-LESIONAL) 1 mm isotropic T1w-like image.
# ------------------------------------------------------------------------------
# OBSERVATION:
# SynthSR improves contrast for surface reconstruction but introduces local
# geometric distortions. Therefore, the synthetic image is non-linearly mapped
# back to the original NLM-denoised anatomy before FastSurfer.

reg_img_sr_seg=${idBIDS}_rec-NLMs2b1_acq-1mm_desc-synthSR_synthseg_UNIT1.nii.gz

mri_synthsr --i ${reg_img_seg_res} --o ${reg_img_sr} --threads 20 --cpu

mri_synthseg --i ${reg_img_sr} --o ${reg_img_sr_seg} --robust --threads 20 --cpu --resample synthseg

# ------------------------------------------------------------------------------
# REGISTRATION 
# ------------------------------------------------------------------------------
# 1. Align SynthSeg resampled image back to native NLM space (rigid).
# 2. Estimate non-linear deformation between SynthSR segmentation and original segmentation.
# 3. Apply transforms to warp SynthSR image into native geometry.
# This preserves SynthSR contrast while recovering original anatomy.

antsRegistrationSyN.sh -d 3 \
    -f ${reg_img} \
    -m ${reg_img_seg_res} \
    -o from-NLMs2b1synth_to-NLMs2b1_ \
    -t r \
    -n 25 \
    -p d


antsRegistrationSyN.sh -d 3 \
    -f ${reg_img_seg} \
    -m ${reg_img_sr_seg} \
    -o from-synthSR_to-NLMs2b1_ \
    -t s \
    -n 25 \
    -p d


antsApplyTransforms -d 3 -u int \
    -i ${reg_img_sr} \
    -r ${reg_img} \
    -t from-NLMs2b1synth_to-NLMs2b1_0GenericAffine.mat \
    -t from-synthSR_to-NLMs2b1_1Warp.nii.gz \
    -t from-synthSR_to-NLMs2b1_0GenericAffine.mat \
    -o sub-PNE032_ses-a1_rec-nlmS2B1_acq-05mm_UNIT1_synthsr.nii.gz


# ------------------------------------------------------------------------------
# 6. FASTSURFER v2.5.4
# ------------------------------------------------------------------------------
#   --vox_size min (Preserve highest available resolution.)
#   --no_biasfield (Bias correction already applied upstream.)
#   --device cpu

t1_2proc=/host/yeatman/local_raid/rcruces/tmp/double_cortex/sub-PNE032_ses-a1_rec-nlmS2B1_acq-05mm_UNIT1_synthsr.nii.gz

dir_surf=/host/yeatman/local_raid/rcruces/tmp/double_cortex/fastsurfer

export SINGULARITY_TMPDIR=/host/bb-compx-03/export02/tmp/singularity_tmp
export SINGULARITY_CACHEDIR=/host/bb-compx-03/export02/tmp/singularity_cache

threads=20

singularity_img=/data_/mica1/01_programs/fastsurfer/fastsurfer-cpu-v2.5.4.sif
fs_licence=/data_/mica1/01_programs/freesurfer-7.3.2/license.txt

singularity run --containall \
    -B ${dir_surf}:${dir_surf} \
    -B ${t1_2proc}:/opt/T1w.nii.gz \
    -B ${fs_licence}:${fs_licence} \
    ${singularity_img} \
        --allow_root \
        --fs_license "${fs_licence}" \
        --t1 /opt/T1w.nii.gz \
        --sid "${idBIDS}" \
        --sd "${dir_surf}" \
        --threads "${threads}" \
        --vox_size min \
        --no_biasfield \
        --no_cereb \
        --no_hypothal \
        --device "cpu"


# ------------------------------------------------------------------------------
# 7. MICAPIPE STRUCTURAL PROCESSING
# ------------------------------------------------------------------------------
#   proc_structural
#   proc_surf
#   post_surf

# ------------------------------------------------------------------------------
# 8. SURFACE-BASED MAPPING
# ------------------------------------------------------------------------------
# A) Get the intracortical surfaces (MPC)
# B) Calculate the laplacian from white matter to the ventricles (exclude subcortical nucleai) 1mm
# C) Mapping: map the qT1 and other quantitative relevant acquisitions the full brain profile (GM-SWM-SWM, deepth) ideally 32k

# ------------------------------------------------------------------------------
# 9. Analisis
# Vertex-wise intensity profile GM/WM segmentation mask
# fMRI mapping
# DM

