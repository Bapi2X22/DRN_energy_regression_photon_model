#!/bin/bash
# run_batch.sh
# FILELIST=$1

# cp -r /eos/user/b/bbapi/Energy_regression/CMSSW_13_3_3/src/test/RecoEgamma-EgammaPhotonProducers/models/ /tmp

# singularity exec /cvmfs/unpacked.cern.ch/registry.hub.docker.com/fastml/triton-torchgeo:22.07-py3-geometric \
# tritonserver --model-repository /tmp/models/ --http-port 9000 --grpc-port 9001 --metrics-port 9002 --allow-http=1 > triton.log 2>&1 &

# cd /eos/user/b/bbapi/Energy_regression/CMSSW_13_3_3/src/condor_scripts/
# export HOME=/afs/cern.ch/user/b/bbapi
# source /cvmfs/cms.cern.ch/cmsset_default.sh
# cmsenv

# for file in $(cat "$FILELIST"); do
#   echo "Running cmsRun on $file"
#   cmsRun Zee_dumper_MINIAOD_MC_cfg_copy_3.py inputFile="$file"
# done

# FILELIST=$1

# # Start triton server
# cp -r /eos/user/b/bbapi/Energy_regression/CMSSW_13_3_3/src/test/RecoEgamma-EgammaPhotonProducers/models/ /tmp

# singularity exec /cvmfs/unpacked.cern.ch/registry.hub.docker.com/fastml/triton-torchgeo:22.07-py3-geometric \
# tritonserver --model-repository /tmp/models/ --http-port 9000 --grpc-port 9001 --metrics-port 9002 --allow-http=1 > triton.log 2>&1 &

# # Setup CMSSW environment
# cd /eos/user/b/bbapi/Energy_regression/CMSSW_13_3_3/src/condor_scripts/
# export HOME=/afs/cern.ch/user/b/bbapi
# source /cvmfs/cms.cern.ch/cmsset_default.sh
# cmsenv

# # Process each file
# for file in $(cat "$FILELIST"); do
#   echo "Running cmsRun on $file"
#   cmsRun Zee_dumper_MINIAOD_MC_cfg_copy_3.py inputFile="$file"

#   # Extract output file name (same as basename of input)
#   outputfile=$(basename "$file")

#   # Check if file exists and is larger than 1MB
#   if [ -f "$outputfile" ]; then
#     filesize=$(stat -c%s "$outputfile")
#     if [ "$filesize" -lt 1000000 ]; then
#       echo "Warning: $outputfile is smaller than 1MB ($filesize bytes). Rerunning cmsRun..."
#       cmsRun Zee_dumper_MINIAOD_MC_cfg_copy_3.py inputFile="$file"
#     fi
#   else
#     echo "Warning: $outputfile not found! Rerunning cmsRun..."
#     cmsRun Zee_dumper_MINIAOD_MC_cfg_copy_3.py inputFile="$file"
#   fi
# done


FILELIST=$1

# Start triton server
cp -r /eos/user/b/bbapi/Energy_regression/CMSSW_13_3_3/src/test/RecoEgamma-EgammaPhotonProducers/models/ /tmp

singularity exec /cvmfs/unpacked.cern.ch/registry.hub.docker.com/fastml/triton-torchgeo:22.07-py3-geometric \
tritonserver --model-repository /tmp/models/ --http-port 9000 --grpc-port 9001 --metrics-port 9002 --allow-http=1 > triton.log 2>&1 &

# Setup CMSSW environment
cd /eos/user/b/bbapi/Energy_regression/CMSSW_13_3_3/src/condor_scripts/
export HOME=/afs/cern.ch/user/b/bbapi
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsenv

# Run cmsRun on each file
for file in $(cat "$FILELIST"); do
  echo "Running cmsRun on $file"
  cmsRun Zee_dumper_MINIAOD_MC_cfg_copy_3.py inputFile="$file"
done

# Check outputs after all runs
for file in $(cat "$FILELIST"); do
  outputfile=$(basename "$file")
  if [ -f "$outputfile" ]; then
    filesize=$(stat -c%s "$outputfile")
    if [ "$filesize" -lt 1000000 ]; then
      echo "Warning: $outputfile is smaller than 1MB ($filesize bytes). Rerunning cmsRun..."
      cmsRun Zee_dumper_MINIAOD_MC_cfg_copy_3.py inputFile="$file"
    fi
  else
    echo "Warning: $outputfile not found! Rerunning cmsRun..."
    cmsRun Zee_dumper_MINIAOD_MC_cfg_copy_3.py inputFile="$file"
  fi
done
