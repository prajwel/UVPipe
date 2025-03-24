#!/bin/bash

# Get the input zip file
zip_file="$1"

# Check if the input file exists and is a zip file
if [[ ! -f "$zip_file" || "${zip_file##*.}" != "zip" ]]; then
    echo "Error: '$zip_file' is not a valid zip file."
    exit 1
fi

# Extract the base name without the extension
base_name=$(basename "$zip_file" .zip)

# Create a temporary directory for extraction
random_string=$(openssl rand -hex 4)
temp_dir="temp_dir_${random_string}"

# Unzip the file into the temporary directory
unzip "$zip_file" -d "$temp_dir" > /dev/null

# Check if unzip was successful
if [ $? -ne 0 ]; then
    echo "Error: Failed to unzip '$zip_file'."
    rm -rf "$temp_dir"
    exit 1
fi

# Create a gzipped tarball with the required naming convention
tarball_name="${base_name}.tar_V9.2.gz"
tar -czf "$tarball_name" -C "$temp_dir" .

# Check if tarball creation was successful
if [ $? -eq 0 ]; then
    echo "Successfully created tarball: $tarball_name"
else
    echo "Error: Failed to create tarball."
    rm -rf "$temp_dir"
    exit 1
fi

# Clean up temporary directory
rm -rf "$temp_dir"

exit 0

