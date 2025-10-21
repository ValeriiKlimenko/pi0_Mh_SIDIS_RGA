#!/bin/bash

# Output file
output_file="nSidis_45_50nA_092025_v2.dat"

# Target directory
target_dir="/cache/clas12/rg-a/production/montecarlo/clasdis_pass2/fa18_inb/"

# Find and list files only, not directories, then write to output file
find "$target_dir" -type f -printf "%f\n" > "$output_file"

echo "File names saved to $output_file"
