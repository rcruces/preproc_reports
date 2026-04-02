# #!/bin/bash

# ------------------------------------------------------------
# MICA-7T: 7t2bids version update workflow
# Last update: March 13, 2025
# ------------------------------------------------------------
#
# Once all the changes in the repository are done:
# -----------------------------------------
# Variables
version="v2.2"
docker_name="micalab/7t2bids"
dcm2bids_img="/data/mica1/01_programs/MICA-7t/7t2bids_${version}.sif"

# Check if the Docker image exists and remove it
if [ -f "$dcm2bids_img" ]; then rm "$dcm2bids_img"; fi

# Change permissions to allow read, write, and execute for all users on all files and directories recursively
chmod aug+rxX -R *
chmod ug+w files/*.*

# Build a Docker image with no cache and tag it as the new version (e.g. ${docker_name}:${version})
docker build --no-cache . -t ${docker_name}:${version}
#docker build . -t ${docker_name}:${version}

# Build a Singularity image from the Docker image
singularity build --fakeroot "${dcm2bids_img}" docker-daemon://${docker_name}:${version}

# Check Softeare versions inside the Singularity container
singularity exec --writable-tmpfs --containall "${dcm2bids_img}" deno --version
singularity exec --writable-tmpfs --containall "${dcm2bids_img}" bc -version
singularity exec --writable-tmpfs --containall "${dcm2bids_img}" jq --version
singularity exec --writable-tmpfs --containall "${dcm2bids_img}" deno run --allow-write -ERN jsr:@bids/validator -V

# -----------------------------------------
# Test singularity image
# -----------------------------------------
# Quick test (7 acquisitions)
#sub="PNA007"
#ses="a2"
# Full test (all acquisitions)
sub="PNE018"
ses="a1"

data_dir="/data_/mica3/BIDS_PNI"
dicoms_dir="${data_dir}/sorted/sub-${sub}_ses-${ses}/dicoms"
bids_dir="${data_dir}/rawdata"
sorted_dir="${data_dir}/sorted/sub-${sub}_ses-${ses}/dicoms_sorted"
#bids_dir="${data_dir}/tmp_test/rawdata"
#sorted_dir="${data_dir}/tmp_test/dicoms_sorted/sub-${sub}_ses-${ses}"

# make sorted directory
if [ ! -d "${sorted_dir}" ]; then mkdir -p "${sorted_dir}"; fi

# PNI: Run the Singularity image
version="v2.2"
dcm2bids_img="/data/mica1/01_programs/MICA-7t/7t2bids_${version}.sif"
singularity run --containall --writable-tmpfs \
    -B "${bids_dir}":"${bids_dir}" \
    -B "${dicoms_dir}":"${dicoms_dir}" \
    -B "${sorted_dir}":"${sorted_dir}" \
    "${dcm2bids_img}" --sub "${sub}" --ses "${ses}" \
    --sorted_dir "${sorted_dir}/sub-${sub}_ses-${ses}" \
    --bids_dir "${bids_dir}" \
    --dicoms_dir "${dicoms_dir}" # --force

# MPN: Run the Singularity image
singularity run --containall --writable-tmpfs \
    -B "${bids_dir}":"${bids_dir}" \
    -B "${sorted_dir}":"${sorted_dir}" \
    "${dcm2bids_img}" --sub "${sub}" --ses "${ses}" \
    --bids_dir "${bids_dir}" \
    --sorted_dir "${sorted_dir}" --dicoms_dir "${sorted_dir}" --force

# -----------------------------------------
# If it works:
# -----------------------------------------

# git add and commit the changes that worked

# Tag the repository with the new version
git tag -f -a ${version}

# Push the tag to the repository
git push origin --tags

# Push Docker image to Docker Hub
docker login --username micalab
docker push ${docker_name}:${version}

# -----------------------------------------
# Make a new release in github repository
# -----------------------------------------
# Print all the commits since the last tag
git log v2.0..HEAD --oneline

# Write the release notes in the release page in github

# GUI https://github.com/MICA-MNI/MICA-7t/tree/main
