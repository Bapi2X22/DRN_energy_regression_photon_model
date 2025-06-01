#!/bin/bash

FILELIST="$1"               # File list passed as first argument
DATASET_NAME="$2"           # Dataset name passed as second argument
output_folder="/eos/user/b/bbapi/Energy_regression/CMSSW_13_3_3/src/condor_scripts/"

if [ ! -f "$FILELIST" ]; then
  echo "Error: File list '$FILELIST' does not exist."
  exit 1
fi

mkdir -p "${output_folder}/${DATASET_NAME}"

# Step 1: Run cmsRun on each input file
for i in $(cat "$FILELIST"); do
  echo "Running cmsRun on $i"
  base_name=$(basename "$i")
  output_file="${output_folder}/${DATASET_NAME}/${base_name}"
  cmsRun Zee_dumper_MINIAOD_MC_cfg_copy_3.py inputFile="$i" datasetname="$DATASET_NAME"
done

# Step 2: Check output file existence and size
for i in $(cat "$FILELIST"); do
  base_name=$(basename "$i")
  output_file="${output_folder}/${DATASET_NAME}/${base_name}"

  if [ -f "$output_file" ]; then
    filesize=$(stat -c%s "$output_file")
    if [ "$filesize" -lt 29360128 ]; then
      echo "Warning: $output_file is smaller than 1MB ($filesize bytes). Rerunning cmsRun..."
      cmsRun Zee_dumper_MINIAOD_MC_cfg_copy_3.py inputFile="$i" datasetname="$DATASET_NAME"
    fi
  else
    echo "Warning: $output_file not found! Rerunning cmsRun..."
    cmsRun Zee_dumper_MINIAOD_MC_cfg_copy_3.py inputFile="$i" datasetname="$DATASET_NAME"
  fi
done

# Step 3: Corruption check (up to 3 attempts)
attempt=1
max_attempts=3
while [ $attempt -le $max_attempts ]; do
  echo "Corruption check attempt $attempt"
  corrupted_files=()

  for i in $(cat "$FILELIST"); do
    base_name=$(basename "$i")
    output_file="${output_folder}/${DATASET_NAME}/${base_name}"

    if [ -f "$output_file" ]; then
      if root -l -b -q "$output_file" 2>&1 | grep -q "probably not closed"; then
        echo "Corrupted file: $output_file"
        corrupted_files+=("$i")
      fi
    fi
  done
