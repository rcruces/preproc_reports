#!/usr/bin/env bash
# Copyright (c) 2025 Raul R. Cruces
# ------------------------------------------------------------
# MICs BIDS Export Script
# ------------------------------------------------------------
#
# Description:
# This script prepares a BIDS-like export of the MICs dataset
# from a raw data directory. It performs the following steps:
#
#   1. Reads a subject list from an Excel file and generates a participants.tsv file compliant with BIDS conventions.
#   2. Copies dataset-level metadata (e.g., CITATION.cff).
#   3. Creates a README file describing the dataset.
#   4. Iterates over participants and processes only ses-01:
#        - Creates the required output directory structure
#        - Symlinks anatomical T1w images
#        - Copies resting-state functional data
#
# The script is designed to be re-runnable and avoids overwriting existing data where possible.
#
# ------------------------------------------------------------
# Requirements:
#   - bash (>= 4.x)
#   - python3
#   - python packages: pandas, openpyxl
#   - rsync
#
# ------------------------------------------------------------
# Inputs:
#   - Raw dataset directory:
#       /data_/mica3/BIDS_MICs/rawdata
#
#   - Excel subject list:
#       ${out}/mica_subjlist_NUS.xlsx
#
#     Expected columns:
#       ID, group, drug_resistant, side, sex,
#       age, onset_age, duration, FBTC, Engel
#
# ------------------------------------------------------------
# Outputs:
#   - BIDS-like directory at:
#       /host/bb-compx-01/export01/data/BIDS_MICs
#
#   - Generated files:
#       participants.tsv
#       README
#       CITATION.cff (copied)
#
# ------------------------------------------------------------
# Notes:
#   - Participant IDs are converted to BIDS format (sub-XXX).
#   - Only session "ses-01" is processed.
#   - Anatomical images (*T1w*) are symlinked to save space.
#   - Functional resting-state files (*task-rest*) are copied.
#   - Missing data are skipped with warnings.
#
# ------------------------------------------------------------

set -euo pipefail

orig="/data_/mica3/BIDS_MICs/rawdata"
out="/host/bb-compx-01/export01/data/BIDS_MICs"
xlsx="${out}/mica_subjlist_NUS.xlsx"
bids=${out}/rawdata

# Timer & Beginning
aloita=$(date +%s)

echo "----------------------------------------"
echo "      Organizing MICs data to share"
echo "----------------------------------------"

# Create output BIDS dir
mkdir -p "${bids}"

# ----------------------------------------
# 1. Create participants.tsv from Excel
# ----------------------------------------
python3 <<EOF
import pandas as pd

xlsx = "${xlsx}"
out_tsv = "${bids}/participants.tsv"

df = pd.read_excel(xlsx)

# Drop column "no" if it exists
if "no" in df.columns:
    df = df.drop(columns=["no"])

# Rename ID -> participant_id
df = df.rename(columns={"ID": "participant_id"})

# Ensure BIDS format (sub- prefix)
df["participant_id"] = df["participant_id"].astype(str).apply(lambda x: f"{x}")

# Save as TSV
df.to_csv(out_tsv, sep="\t", index=False)
EOF

echo "Created participants.tsv"

# Create participants.json
# Create participants.json
cat <<EOF > "${bids}/participants.json"
{
  "participant_id": {
    "Description": "Unique participant identifier in BIDS format (sub-XXX)."
  },
  "group": {
    "Description": "Participant group classification.",
    "Levels": {
      "control": "Healthy control participant",
      "TLE": "Temporal lobe epilepsy patient"
    }
  },
  "drug_resistant": {
    "Description": "Indicates whether the participant has drug-resistant epilepsy.",
    "Levels": {
      "0": "No",
      "1": "Yes"
    }
  },
  "side": {
    "Description": "Lateralization of the epileptic focus.",
    "Levels": {
      "L": "Left hemisphere",
      "R": "Right hemisphere",
      "B": "Bilateral",
      "U": "Unknown or not specified",
      "C": "Control (no lateralization)"
    }
  },
  "sex": {
    "Description": "Biological sex of the participant.",
    "Levels": {
      "M": "Male",
      "F": "Female"
    }
  },
  "age": {
    "Description": "Age of the participant at the time of scanning.",
    "Units": "years"
  },
  "onset_age": {
    "Description": "Age at epilepsy onset.",
    "Units": "years"
  },
  "duration": {
    "Description": "Duration of epilepsy at the time of scanning (age minus onset_age).",
    "Units": "years"
  }
}
EOF

# ----------------------------------------
# 2. Copy CITATION.cff
# ----------------------------------------
if [ -f "${orig}/CITATION.cff" ]; then
    cp "${orig}/CITATION.cff" "${bids}/"
    echo "Copied CITATION.cff"
else
    echo "Warning: CITATION.cff not found in ${orig}"
fi

# ----------------------------------------
# 3. Create README
# ----------------------------------------
cat <<EOF > "${bids}/README"
The Microstructure-Informed Connectomics (MICs) dataset was provided by the Multimodal Imaging and Connectome Analysis Lab, Montreal Neurological Institute and Hospital, The Neuro.

If you reference this dataset in your publications, please acknowledge its authors.
EOF

echo "Created README"

# ----------------------------------------
# 4. Loop over participants (ses-01 only)
# ----------------------------------------
participants=$(cut -f1 "${bids}/participants.tsv" | tail -n +2)

for sub in ${participants}; do
    ses="01"

    dir_orig="${orig}/${sub}/ses-${ses}"
    dir_out="${bids}/${sub}/ses-${ses}"

    echo -e "\n>> Processing ${sub} ses-${ses}"

    mkdir -p "${dir_out}/anat" "${dir_out}/func"

    # ------------------------------------
    # Symlink T1w
    # ------------------------------------
    if compgen -G "${dir_orig}/anat/*T1w*" > /dev/null; then
        for f in ${dir_orig}/anat/*T1w*; do
            ln -sf "${f}" "${dir_out}/anat/"
        done
    else
        echo "No T1w found for ${sub}"
    fi

    # ------------------------------------
    # Copy resting-state func
    # ------------------------------------
    if compgen -G "${dir_orig}/func/*task-rest*" > /dev/null; then
        rsync -av --ignore-existing \
            ${dir_orig}/func/*task-rest* \
            "${dir_out}/func/"
    else
        echo "No rest func found for ${sub}"
    fi

done

# -----------------------------------------------------------------------------------------------#
# Rename sbref AP PA to be bids compliant

# ANONIMIZATION of T1w data

# -----------------------------------------------------------------------------------------------#
# BIDS validation
echo "bids_validator_output.txt" >> ${bids}/.bidsignore
deno run --allow-write -ERN jsr:@bids/validator ${bids} --ignoreWarnings --outfile ${bids}/bids_validator_output.txt

# -----------------------------------------------------------------------------------------------#
#			 Total Running Time
lopuu=$(date +%s)
eri=$(echo "$lopuu - $aloita" | bc)
eri=$(echo print $eri/60 | perl)
Title "TOTAL running time:\033[38;5;220m $(printf "%0.3f\n" ${eri}) minutes \033[38;5;141m"