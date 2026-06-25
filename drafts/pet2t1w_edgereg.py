#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Edge-based affine PET → T1w registration.
 
Steps:
  0. Load inputs          ants.image_read (nibabel for output header only)
  1. Smooth both          ants.smooth_image (sigma=1 mm physical)
  2. Edge filter both     np.gradient magnitude scaled by voxel spacing
  3. Affine registration  ants.registration with internal histogram matching (-j 1)
  4. Apply transform      ants.apply_transforms on original PET

Created on June 2026, modified from bash workflow by rcruces
@author: rcruces
"""

import numpy as np
import nibabel as nib
import ants
 
 
def pet_to_t1w_edge_registration(
    pet_path: str,
    t1w_path: str,
    output_prefix: str,
    smooth_sigma_mm: float = 1.0,
    interpolator: str = "nearestNeighbor",
) -> nib.Nifti1Image:
    """
    Register PET to T1w using edge-based affine registration.
 
    Mirrors: antsRegistrationSyN.sh -d 3 -t a -j 1 on smoothed edge images,
    then antsApplyTransforms -n NearestNeighbor on the original PET.
 
    Parameters
    ----------
    pet_path        : path to PET NIfTI (average or single frame)
    t1w_path        : path to T1w NIfTI (registration target)
    output_prefix   : ANTs output prefix; affine .mat saved here
    smooth_sigma_mm : Gaussian sigma in mm (default 1.0 → fslmaths -s 1)
    interpolator    : interpolation for final warp (default "nearestNeighbor")
 
    Returns
    -------
    nib.Nifti1Image  PET resampled into T1w space
    """
 
    # ------------------------------------------------------------------
    # Step 0 – load
    t1w_nib = nib.load(t1w_path)          # kept only for output header/affine
    pet = ants.image_read(pet_path)
    t1w = ants.image_read(t1w_path)
 
    # ------------------------------------------------------------------
    # Step 1 – smooth  (fslmaths -s sigma)
    pet_smooth = ants.smooth_image(pet, smooth_sigma_mm, sigma_in_physical_coordinates=True)
    t1w_smooth = ants.smooth_image(t1w, smooth_sigma_mm, sigma_in_physical_coordinates=True)
 
    # ------------------------------------------------------------------
    # Step 2 – edge filter  (fslmaths -edge → gradient magnitude)
    def _edge(img):
        d, sp = img.numpy().astype(np.float64), np.array(img.spacing)
        g = [np.gradient(d, sp[i], axis=i) for i in range(3)]
        return img.new_image_like(np.sqrt(sum(gi**2 for gi in g)))
 
    pet_edge = _edge(pet_smooth)
    t1w_edge = _edge(t1w_smooth)
 
    # ------------------------------------------------------------------
    # Step 3 – affine registration on edge images
    # use_legacy_histogram_matching=True mirrors -j 1 in antsRegistrationSyN.sh
    reg = ants.registration(
        fixed             = t1w_edge,
        moving            = pet_edge,
        type_of_transform = "Affine",
        outprefix         = output_prefix,
        use_legacy_histogram_matching = True,
    )
 
    # ------------------------------------------------------------------
    # Step 4 – apply transform to original (unmodified) PET
    pet_warped = ants.apply_transforms(
        fixed         = t1w,
        moving        = pet,
        transformlist = reg["fwdtransforms"],
        interpolator  = interpolator,
    )
 
    return nib.Nifti1Image(pet_warped.numpy(), t1w_nib.affine, t1w_nib.header)

# ----------------------------------------------------------------------
# CLI
if __name__ == "__main__":
    import argparse
 
    parser = argparse.ArgumentParser(description="Edge-based affine PET to T1w registration")
    parser.add_argument("pet", help="path to PET NIfTI file")
    parser.add_argument("t1w", help="path to T1w NIfTI file")
    parser.add_argument("output_prefix", help="prefix for ANTs output files (affine .mat saved here)")
    parser.add_argument("--sigma", type=float, default=1.0, help="smoothing sigma in mm (default: 1.0)")
    parser.add_argument("--interpolator", default="nearestNeighbor", help="interpolation for final warp (default: nearestNeighbor)")
    args = parser.parse_args()
 
    result = pet_to_t1w_edge_registration(
        args.pet, args.t1w, args.output_prefix,
        smooth_sigma_mm=args.sigma,
        interpolator=args.interpolator,
    )
    out = args.output_prefix + "_edgereg.nii.gz"
    nib.save(result, out)
    print(f"Saved: {out}")

