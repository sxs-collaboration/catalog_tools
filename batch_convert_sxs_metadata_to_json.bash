#!/bin/bash

# Distributed under the MIT License.
# See LICENSE.txt for details.

# This script takes as arguments one or more paths. Each path provided
# should contain a file called metadata.txt, and each path should end in 
# "/Lev#" or "/Res#", where # is an integer. For each path, the script will 
# convert the metadata.txt file to metadata.json.
#
# To use this script, first edit the following parameters

# Path to conversion script
convert_path=/home/geoffrey/BBH/Catalog/catalog_tools/convert_sxs_metadata_txt_to_json.py

# Batch command used when submitting jobs to convert waveforms
batch_cmd="srun -p orca-0 -n 1 -t 12:00:00"

# Python environment command
# Run this command to activate the python environment you want to use
# Set to "echo date" or some other do-nothing command if this isn't 
# necessary on your system.
python_env_cmd="/home/geoffrey/apps/anaconda3/bin/activate root"

#############################################################
# SHOULD NOT NEED TO MODIFY ANYTHING BELOW THIS LINE
# (except, possibly, flags in cmd_to_run)
#############################################################

# First, make sure we can find the conversion script
convert_cmd=$(readlink -f ${convert_path})
if [[ -f ${convert_cmd} ]]
then
    echo "Found conversion script ${convert_cmd}"
else
    echo "ERROR: could not find ${convert_cmd}"
    exit 1
fi

# Check that conversion command works
source ${python_env_cmd}
echo "Using $(which python)"
echo "Testing python ${convert_cmd} --help"
$(python ${convert_cmd} --help &> /dev/null)
convert_test_status=$?
if [ ${convert_test_status} -ne 0 ]
then
    echo "ERROR: could not run conversion script"
fi

# Next, check that all paths have the expected data and name
for path in "$@"
do
    abs_path=$(readlink -f ${path})
    echo "Checking ${abs_path}"

    # Get resolution and paths to data to read
    base_path=$(basename ${path})
    resolution=$(echo ${base_path} | sed "s/Lev//" | sed "s/Res//")
    metadata_file=${abs_path}/metadata.txt

    if [ "${resolution}" -eq "${resolution}" ]
    then
        echo "Resolution: ${resolution} (inferred from ${base_path})"
    else
        echo "ERROR: cannot infer resolution from path ${base_path}"
        exit 1
    fi

    # Check that required data files exist
    files_to_check=(${metadata_file})
    for file in ${files_to_check[@]}
    do
        if [[ -f "${file}" ]]
        then
            echo "Found $file"
        else
            echo "ERROR: Could not find ${file}"
            exit 1
        fi
    done
done

# Finally, loop over all paths and submit a conversion job
for path in "$@"
do
    abs_path=$(readlink -f ${path})
    echo "Processing ${abs_path}"
    base_path=$(basename ${path})
    resolution=$(echo ${base_path} | sed "s/Lev//" | sed "s/Res//")
    metadata_file=${abs_path}/metadata.txt
    
    echo "Running: python ${convert_cmd} --txt ${metadata_file} --output ${abs_path}/metadata.json"
    python ${convert_cmd} --txt ${metadata_file} --output ${abs_path}/metadata.json
done
