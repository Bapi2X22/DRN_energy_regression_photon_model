#!/bin/bash

# Loop over runlist indices
for i in {1..4}; do
    # Generate the runlist file name
    runlist="runlist_${i}.txt"
        
    # Generate the output file name
    output_file="Plot_final_barrel_${i}.root"
    
    # Run the command
    ./analyzeHGCMuons_D1 "$runlist" "$output_file" data 2500
done
