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
        return f"{directory_prefix}/" + "_join(filename_parts)" + f"_{suffix}"

class BIDS_micr_metadata:
    """Extracts and stores BIDS metadata for microscopy."""
    def __init__(self, ndpi_path):
        self.path = ndpi_path
        self.metadata = {
            "PixelSize": None, "PixelSizeUnits": "um", "Immersion": "Air",
            "NumericalAperture": None, "Magnification": None,
            "ImageAcquisitionProtocol": "Antigen retrieval CC1 buffer (24 min at 91 °C) on Ventana Discovery Ultra",
            "OtherAcquisitionParameters": "Digitized at 40x; interpreted by expert neuropathologist",
            "BodyPart": "BRAIN", 
            "BodyPartDetails": "R/L cortex/hippocampus/amigdala", 
            "BodyPartDetailsOntology": "https://www.ebi.ac.uk/ols/ontologies/uberon",
            "SampleEnvironment": "ex vivo", "SampleEmbedding": "paraffin", "SampleFixation": "formalin 4%",
            "SampleStaining": "Immunohistochemistry for phosphorylated tau", 
            "SamplePrimaryAntibody": "AT8 (Ser202/Thr205; Thermo Fisher Scientific; MN1020); dilution 1:1200", 
            "Manufacturer": "Hamamatsu", "ManufacturersModelName": "NanoZoomer S210", 
            "InstitutionName": "Montreal Neurological Institute, McGill University", 
            "InstitutionAddress": "3801 University St, Montreal, Quebec H3A 2B4, Canada", 
            "InstitutionalDepartmentName": "Neuropathology"
        }

    def fill_from_ndpi(self):
        with tifffile.TiffFile(self.path) as tif:
            if tif.ome_metadata:
                root = ET.fromstring(tif.ome_metadata)
                ns = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06'}
                pixels = root.find('.//ome:Pixels', ns)
                if pixels is not None:
                    self.metadata["PixelSize"] = [
                        float(pixels.get('PhysicalSizeX', 0)),
                        float(pixels.get('PhysicalSizeY', 0))
                    ]
                obj = root.find('.//ome:Objective', ns)
                if obj is not None:
                    self.metadata["Magnification"] = float(obj.get('NominalMagnification', 40))
                    self.metadata["NumericalAperture"] = float(obj.get('LensNA', 0))

            ndpi_info = getattr(tif, 'ndpi_tags', {})
            self.metadata["DeviceSerialNumber"] = ndpi_info.get('SerialNumber', 'Unknown')
            self.metadata["SoftwareVersions"] = ndpi_info.get('SoftwareVersion', 'Unknown')

    def save_json(self, output_path):
        with open(output_path, 'w') as f:
            json.dump(self.metadata, f, indent=4)

# --- Execution Script ---

def main():
    parser = argparse.ArgumentParser(description="Convert NDPI to BIDS Microscopy structure.")
    
    # Mandatory Arguments
    parser.add_argument("--ndpi_path", required=True, help="Path to source .ndpi file")
    parser.add_argument("--bids_root", required=True, help="Path to BIDS dataset root")
    parser.add_argument("--sub", required=True, help="Subject ID (e.g., 01)")
    parser.add_argument("--stain", required=True, help="Stain label (e.g., AT8)")
    parser.add_argument("--suffix", required=True, help="BIDS suffix (e.g., BF)")

    # Optional BIDS Entities
    parser.add_argument("--ses", help="Session ID")
    parser.add_argument("--sample", default="01", help="Sample ID")
    parser.add_argument("--acq", help="Acquisition label (e.g., 40x)")
    parser.add_argument("--run", help="Run index")
    parser.add_argument("--chunk", help="Chunk index")
    parser.add_argument("--convert", action="store_true", help="Convert to OME-TIFF using bfconvert")

    args = parser.parse_args()

    # 0. Setup Logging
    log_file = os.path.join(args.bids_root, "bids_conversion.log")
    os.makedirs(args.bids_root, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()]
    )

    logging.info(f"Processing started for: {args.ndpi_path}")

    # 1. Setup Naming Entities
    entities = vars(args)
    bids_namer = BIDS_micr_name(**entities)
    bids_rel_path = bids_namer.build()
    full_bids_base = os.path.join(args.bids_root, bids_rel_path)

    try:
        # 2. Create Directory
        os.makedirs(os.path.dirname(full_bids_base), exist_ok=True)

        # 3. Copy NDPI
        target_ndpi = f"{full_bids_base}.ndpi"
        logging.info(f"Copying NDPI to {target_ndpi}")
        shutil.copy2(args.ndpi_path, target_ndpi)

        # 4. Metadata Extraction
        logging.info("Extracting metadata...")
        meta = BIDS_micr_metadata(target_ndpi)
        meta.fill_from_ndpi()
        meta.save_json(f"{full_bids_base}.json")
        logging.info(f"Metadata JSON saved.")

        # 5. Optional OME-TIFF Conversion
        if args.convert:
            target_ome = f"{full_bids_base}.ome.tif"
            logging.info("Starting OME-TIFF conversion...")
            result = subprocess.run(
                ['bfconvert', '-bigtiff', '-compression', 'LZW', target_ndpi, target_ome], 
                capture_output=True, text=True
            )
            if result.returncode == 0:
                logging.info(f"Conversion successful: {target_ome}")
            else:
                logging.error(f"bfconvert failed: {result.stderr}")

        logging.info(f"Successfully processed sub-{args.sub}")

    except Exception as e:
        logging.error(f"Failed to process {args.ndpi_path}: {str(e)}")

if __name__ == "__main__":
    main()