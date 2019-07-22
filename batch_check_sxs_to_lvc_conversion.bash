#!/bin/bash
#!/bin/bash -
#SBATCH -o check_sxs_to_lvc.stdout
#SBATCH -e check_sxs_to_lvc.stderr
#SBATCH --ntasks 20
#SBATCH --cpus-per-task 1
#SBATCH -J check_sxs_to_lvc
#SBATCH --nodes 1
#SBATCH -p orca-1
#SBATCH -t 48:00:00
#SBATCH -D /home/geoffrey/BBH/Catalog/convert_2018_public

# Distributed under the MIT License.
# See LICENSE.txt for details.

# This script checks the conversion of an SXSsimulation to LVC format.
# Checks include the following:
#    1. Run lvcnrcheck on the LVC-format file
#    2. Run compare_sxs_lvc.py, to check that
#       the LVC data agrees with the SXS data
#
# The script will check all runs found in path_to_check for directories
# named ${level_dirname}*
path_to_check="/home/geoffrey/BBH/Catalog/zenodo"
level_dirname="Lev"

# Python environment command
# Run this command to activate the python environment you want to use
# Set to "echo date" or some other do-nothing command if this isn't 
# necessary on your system.
python_env_cmd="/home/geoffrey/apps/anaconda2/bin/activate root"

# lvcnrcheck path
# This should point to a copy of lvcnrcheck, which is part of lvcnrpy
# https://git.ligo.org/waveforms/lvcnrpy
lvcnrcheck_path="/home/geoffrey/BBH/lvcnrpy/bin/lvcnrcheck"

# LVCNR_DATADIR
# This shell variable should point to a checkout of the ligo lvc git repo
export LVCNR_DATADIR=/home/geoffrey/BBH/lvcnr-lfs

# Path to compare_sxs_vs_lvc.py, which is part of 
# https://github.com/sxs-collaboration/catalog_tools
compare_sxs_vs_lvc_path="/home/geoffrey/BBH/Catalog/catalog_tools/compare_sxs_vs_lvc.py"

# Romspline path...should point to a directory named romspline, containing
# a checkout of
# https://bitbucket.org/chadgalley/romspline/src/master/
# (e.g., this is included in private SXS repo 
# https://github.com/sxs-collaboration/CatalogAnalysis)
romspline_path="/home/geoffrey/BBH/Catalog/CatalogAnalysis"

#############################################################
# SHOULD NOT NEED TO MODIFY ANYTHING BELOW THIS LINE
#############################################################
# Find  all the directories to check
levels=$(find ${path_to_check} -name ${level_dirname}*|xargs)
for level in $levels
do
    echo "Found directory to check: $level"
done

# Make sure we can find the check scripts
lvcnrcheck_cmd=$(readlink -f ${lvcnrcheck_path})
if [[ -f ${lvcnrcheck_cmd} ]]
then
    echo "Found lvcnrcheck script ${lvcnrcheck_cmd}"
else
    echo "ERROR: could not find ${lvcnrcheck_cmd}"
    exit 1
fi

compare_sxs_vs_lvc_cmd=$(readlink -f ${compare_sxs_vs_lvc_path})
if [[ -f ${compare_sxs_vs_lvc_cmd} ]]
then
    echo "Found compare_sxs_vs_lvc script ${compare_sxs_vs_lvc_cmd}"
else
    echo "ERROR: could not find ${compare_sxs_vs_lvc_cmd}"
    exit 1
fi

# Check that check commands work
source ${python_env_cmd}
echo "Using $(which python)"
echo "Testing python ${compare_sxs_vs_lvc_cmd} --help"
$(python ${compare_sxs_vs_lvc_cmd} --help &> /dev/null)
convert_test_status=$?
if [ ${convert_test_status} -ne 0 ]
then
    echo "ERROR: could not run compare_sxs_vs_lvc.py"
fi

echo "Testing lvcnrcheck --help"
$(python ${lvcnrcheck_cmd} --help &> /dev/null)
convert_test_status=$?
if [ ${convert_test_status} -ne 0 ]
then
    echo "ERROR: could not run lvcnrcheck"
fi


# Next, check that all paths have the expected data and name
for path in $levels
do
    abs_path=$(readlink -f ${path})
    echo "Checking ${abs_path}"

    # Get resolution and paths to data to read
    base_path=$(basename ${path})
    resolution=$(echo ${base_path} | sed "s/Lev//" | sed "s/Res//")
    waveform_file=${abs_path}/rhOverM_Asymptotic_GeometricUnits_CoM.h5
    horizons_file=${abs_path}/Horizons.h5
    metadata_file=${abs_path}/metadata.json
    lvc_filename=$(basename $(dirname ${abs_path}))_Res${resolution}.h5
    lvc_file=${abs_path}/${lvc_filename}

    if [ "${resolution}" -eq "${resolution}" ]
    then
        echo "Resolution: ${resolution} (inferred from ${base_path})"
    else
        echo "ERROR: cannot infer resolution from path ${base_path}"
        exit 1
    fi

    # Check that required data files exist
    files_to_check=(${waveform_file} ${horizons_file} ${metadata_file} ${lvc_file})
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

for path in $levels
do
    abs_path=$(readlink -f ${path})
    echo "Processing ${abs_path}"

    # Get resolution and paths to data to read
    base_path=$(basename ${path})
    resolution=$(echo ${base_path} | sed "s/Lev//" | sed "s/Res//")
    waveform_file=${abs_path}/rhOverM_Asymptotic_GeometricUnits_CoM.h5
    horizons_file=${abs_path}/Horizons.h5
    metadata_file=${abs_path}/metadata.json
    lvc_filename=$(basename $(dirname ${abs_path}))_Res${resolution}.h5
    lvc_file=${abs_path}/${lvc_filename}

    echo "Running lvcnrcheck on ${abs_path}"
    ${lvcnrcheck_cmd} -f 3 ${lvc_file}

    echo "Running compare_sxs_vs_lvc.py on ${abs_path}"
    python ${compare_sxs_vs_lvc_cmd} --lvc_file ${lvc_file} \
                                     --sxs_waveform ${waveform_file} \
                                     --sxs_horizons ${horizons_file} \
                                     --sxs_json ${metadata_file} \
                                     --romspline_path ${romspline_path}
done
