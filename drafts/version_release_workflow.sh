# #!/bin/bash

# ------------------------------------------------------------
# MICA-7T: 7t2bids version update workflow
# Last update: March 13, 2025
# ------------------------------------------------------------

# Once all the changes in the repository are done:

# Tag the repository with the new version
git tag -a v2.2 -m "Version 2.2"

# Push the tag to the repository
git push origin --tags

# Path to Singularity
dcm2bids_img="/data/mica1/01_programs/MICA-7t/7t2bids_v2.2.sif"

# Change permissions to allow read, write, and execute for all users on all files and directories recursively
chmod aug+rxX -R *

# Build a Docker image with no cache and tag it as 7t2bids:v2.2
docker build --no-cache . -t 7t2bids:v2.2

# Build a Singularity image from the Docker image
singularity build "${dcm2bids_img}" docker-daemon://7t2bids:v2.2

# Check the version of deno inside the Singularity container
singularity exec --writable-tmpfs --containall "${dcm2bids_img}" deno --version

# Check the version of bc inside the Singularity container
singularity exec --writable-tmpfs --containall "${dcm2bids_img}" bc -version

# Check the version of jq inside the Singularity container
singularity exec --writable-tmpfs --containall "${dcm2bids_img}" jq --version

# Check the version of bids validator inside the Singularity container
singularity exec --writable-tmpfs --containall "${dcm2bids_img}" deno run --allow-write -ERN jsr:@bids/validator -V

# ------------------------------
# Test singularity image
# ------------------------------
# quick test (7 acquisitions)
# sub="PNA007"
# ses="a2"

# full test (all acquisitions)
sub="PNC019"
ses="a1"
data_dir="/data_/mica3/BIDS_PNI"
bids_dir="${data_dir}/tmp_test/rawdata"
dicoms_dir="${data_dir}/sorted/sub-${sub}_ses-${ses}/dicoms"
sorted_dir="${data_dir}/tmp_test/dicoms_sorted/sub-${sub}_ses-${ses}"

# make sorted directory
mkdir -p "${sorted_dir}"

# Run the Singularity image
singularity run --containall \
    -B "${bids_dir}":"${bids_dir}" \
    -B "${dicoms_dir}":"${dicoms_dir}" \
    -B "${sorted_dir}":"${sorted_dir}" \
    "${dcm2bids_img}" --sub "${sub}" --ses "a3" \
    --dicoms_dir "${dicoms_dir}" \
    --sorted_dir "${sorted_dir}/sub-${sub}_ses-a3" \
    --bids_dir "${bids_dir}"

# If it works:
# Push Docker image to Docker Hub
docker tag 7t2bids:v2.2 mica7t/7t2bids:v2.2
docker login --username micalab
docker push mica7t/7t2bids:v2.2

# ------------------------------
# Make a new release in github repository
# ------------------------------
# Print all the commits since the last tag
git log v2.0..HEAD --oneline

# Write the release notes in the release page in github

# GUI https://github.com/MICA-MNI/MICA-7t/tree/main
