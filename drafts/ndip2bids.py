#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on March 12, 2026
@author: rcruces
"""

import os
import json
import shutil
import argparse
import subprocess
import tifffile
import logging
import xml.etree.ElementTree as ET

# --- Class Definitions ---

class BIDS_micr_name:
    """Generates BIDS-compliant paths and filenames for microscopy data."""
    def __init__(self, **kwargs):
        self.entities = ["sub", "ses", "sample", "acq", "stain", "run", "chunk"]
        self.values = kwargs

    def build(self):
        sub_val = self.values.get("sub")
        ses_val = self.values.get("ses")
        
        path_segments = [f"sub-{sub_val}"]
        if ses_val:
            path_segments.append(f"ses-{ses_val}")
        path_segments.append("micr")
        
        directory_prefix = "/".join(path_segments)

        filename_parts = []
        for entity in self.entities:
            if entity in self.values and self.values[entity]:
                filename_parts.append(f"{entity}-{self.values[entity]}")

        suffix = self.values.get("suffix", "BF")
        filename = "_".join(filename_parts) + f"_{suffix}"
        return f"{directory_prefix}/{filename}"

class BIDS_micr_metadata:
    """Extracts and stores BIDS metadata for microscopy."""
    def __init__(self, ndpi_path):
        self.path = ndpi_path
        self.metadata = {            
            # Get from metadata
            # Note pixel units MUST be µm
            "Manufacturer": None,
            "ManufacturersModelName": None,
            "PixelSize": None,
            "PixelSizeUnits": None,
            "Magnification": None,
            "ImageAcquisitionProtocol": None,
            "ScanTimeSeconds": None,
            "FocusTimeSeconds": None,
            "Software": None,
            "DateAcquired": None,
            "Compression": None,
            "BitsPerPixel": None,

            # Institution
            "InstitutionName": "Montreal Neurological Institute, McGill University", 
            "InstitutionAddress": "3801 University St, Montreal, Quebec H3A 2B4, Canada", 
            "InstitutionalDepartmentName": "Neuropathology"
        }

    def fill_from_ndpi(self):
        """Deep search extraction for stubborn metadata fields."""
        with tifffile.TiffFile(self.path) as tif:
            page = tif.pages[0]
            
            # 1. Standard OME-XML Wildcard
            try:
                if tif.ome_metadata:
                    root = ET.fromstring(tif.ome_metadata.strip())
                    pixels = root.find(".//{*}Pixels")
                    if pixels is not None:
                        self.metadata["PixelSize"] = [
                            round(float(pixels.get('PhysicalSizeX', 0)), 4),
                            round(float(pixels.get('PhysicalSizeY', 0)), 4)
                        ]
                    
                    obj = root.find(".//{*}Objective")
                    if obj is not None:
                        na = obj.get('LensNA')
                        if na: self.metadata["NumericalAperture"] = float(na)
            except: pass

            # 2. Proprietary NDPI Tags Search
            ndpi_info = getattr(tif, 'ndpi_tags', {})
            if ndpi_info:
                # Some scanners store it directly in this dict
                if not self.metadata["NumericalAperture"]:
                    self.metadata["NumericalAperture"] = ndpi_info.get('NA') or ndpi_info.get('NumericalAperture')
                
                if not self.metadata["PixelSize"] and 'Distance' in ndpi_info:
                    psize = round(float(ndpi_info['Distance']) / 1000.0, 4)
                    self.metadata["PixelSize"] = [psize, psize]

            # 3. Broader Regex for ImageDescription
            # Handles "NA: 0.75", "N.A. 0.75", "NumericalAperture=0.75", etc.
            desc = page.tags.get('ImageDescription')
            if desc and not self.metadata["NumericalAperture"]:
                desc_str = str(desc.value)
                # Regex: Look for variations of NA followed by a float
                na_match = re.search(r'(?:NA|N\.A\.|Numerical\s?Aperture)[:\s=]+([0-9\.]+)', desc_str, re.IGNORECASE)
                if na_match:
                    self.metadata["NumericalAperture"] = float(na_match.group(1))

            # 4. Brute-force TIFF Tag Scan (Last Resort)
            if not self.metadata["NumericalAperture"]:
                for tag in page.tags:
                    tag_data = str(tag.value)
                    if "NA" in tag_data or "Aperture" in tag_data:
                        na_match = re.search(r'([0-9]\.[0-9]{2})', tag_data)
                        if na_match:
                            self.metadata["NumericalAperture"] = float(na_match.group(1))
                            break

            # 5. Resolution Fallback for PixelSize
            if not self.metadata["PixelSize"]:
                x_res = page.tags.get('XResolution')
                unit = page.tags.get('ResolutionUnit')
                if x_res and unit:
                    res_val = x_res.value[0] / x_res.value[1]
                    if unit.value == 3: psize = 10000.0 / res_val
                    elif unit.value == 2: psize = 25400.0 / res_val
                    self.metadata["PixelSize"] = [round(psize, 4), round(psize, 4)]

    def save_json(self, output_path):
        """Writes the dictionary to a BIDS-compliant JSON file."""
        with open(output_path, 'w') as f:
            json.dump(self.metadata, f, indent=4)

# --- Execution Script ---

def main():
    parser = argparse.ArgumentParser(description="Convert NDPI files to BIDS Microscopy format (MNI Neuropathology).")
    
    # Mandatory BIDS arguments
    parser.add_argument("--ndpi_path", required=True, help="Path to raw Hamamatsu .ndpi file")
    parser.add_argument("--bids", required=True, help="Path to the root of the BIDS dataset")
    parser.add_argument("--sub", required=True, help="Subject ID (e.g., PX067)")
    parser.add_argument("--stain", required=True, help="Stain entity (e.g., AT8)")
    parser.add_argument("--suffix", required=True, help="BIDS suffix (e.g., BF)")

    # Optional BIDS entities
    parser.add_argument("--ses", help="Session ID")
    parser.add_argument("--sample", help="Sample ID (e.g., NP24709)")
    parser.add_argument("--acq", help="Acquisition label")
    parser.add_argument("--run", help="Run index")
    parser.add_argument("--chunk", help="Chunk label (e.g., A3)")
    
    # Operational flags
    parser.add_argument("--convert", action="store_true", help="Convert NDPI to OME-TIFF via bfconvert")
    parser.add_argument("--force", action="store_true", help="Overwrite existing BIDS files and sidecars")

    args = parser.parse_args()

    # --- 1. Setup Naming and Directory Logic ---
    entities = {k: v for k, v in vars(args).items() if v is not None}
    bids_namer = BIDS_micr_name(**entities)
    bids_rel_path = bids_namer.build()
    full_bids_base = os.path.join(args.bids, bids_rel_path)
    
    # Define log location at the same level as /micr
    subject_session_root = os.path.dirname(os.path.dirname(full_bids_base))
    log_dir = os.path.join(subject_session_root, "log")
    os.makedirs(log_dir, exist_ok=True)

    # --- 2. Setup Logging ---
    log_name_parts = [f"sub-{args.sub}"]
    if args.ses: log_name_parts.append(f"ses-{args.ses}")
    log_name_parts.append("bids-micr.log")
    
    log_file = os.path.join(log_dir, "_".join(log_name_parts))
    
    # Configure logging to write to the new subject-specific file
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()]
    )

    logging.info(f"Processing started for: {args.ndpi_path}")

    try:
        # Step 3: Directories and Copy
        os.makedirs(os.path.dirname(full_bids_base), exist_ok=True)
        target_ndpi = f"{full_bids_base}.ndpi"
        target_json = f"{full_bids_base}.json"

        if os.path.exists(target_ndpi) and not args.force:
            logging.info(f"SKIPPED: {target_ndpi} exists. Use --force to overwrite.")
            return

        status_prefix = "OVERWRITE" if os.path.exists(target_ndpi) else "NEW"
        logging.info(f"{status_prefix}: Copying NDPI to {target_ndpi}")
        shutil.copy2(args.ndpi_path, target_ndpi)

        # Step 4: Metadata
        logging.info("Metadata: Extracting from headers...")
        meta = BIDS_micr_metadata(target_ndpi)
        meta.fill_from_ndpi()
        meta.save_json(target_json)
        logging.info("Metadata: JSON sidecar saved.")

        # Step 5: Optional Conversion
        if args.convert:
            target_ome = f"{full_bids_base}.ome.tif"
            if os.path.exists(target_ome) and not args.force:
                logging.info("Conversion: OME-TIFF exists, skipping.")
            else:
                logging.info("Conversion: Starting bfconvert...")
                result = subprocess.run(
                    ['bfconvert', '-bigtiff', '-compression', 'LZW', target_ndpi, target_ome],
                    capture_output=True, text=True
                )
                if result.returncode == 0:
                    logging.info(f"Conversion: Success -> {target_ome}")
                else:
                    logging.error(f"Conversion: Failed -> {result.stderr}")

        logging.info(f"STATUS: SUCCESS for sub-{args.sub}")

    except Exception as e:
        logging.error(f"STATUS: FAILED for sub-{args.sub} - {str(e)}")

if __name__ == "__main__":
    main()