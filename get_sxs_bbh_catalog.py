# Distributed under the MIT License.
# See LICENSE.txt for details.

# A script that downloads the following for all resolutions
# of all binary-black-hole simulations in the publc SXS catalog,
# placing the results in the directory specified by the 
# --output_dir option:
#
#  * metadata.txt (converted to metadata.json)
#  * rhOverM_Asymptotic_GeometricUnits_CoM.h5
#  * Horizons.h5

import argparse
import json
import os
import sxs
from sxs import zenodo as zen


def bbh_keys_from_simulation_keys(simulation_keys):
    bbh_simulations = []
    for simulation in simulation_keys:
        if simulation.split(':')[-2] == "BBH":
            bbh_simulations.append(simulation)
    return bbh_simulations

p = argparse.ArgumentParser(description="Get SXS public BBH data")
p.add_argument("--output_dir",
               help="Name of output directory (must not exist)",
               required=True)
p.add_argument("--sxs_catalog_metadata",
               help="Path to sxs_catalog.json (via get_sxs_public_metadata.py)",
               required=True)
p.add_argument("--highest_lev_only",
               help="Only download data for the highest available resolution",
               action="store_true")
p.add_argument("--dry_run",
               help="Test this script, but don't actually download anything",
               action="store_true")
args = p.parse_args()

try:
    os.mkdir(args.output_dir)
except:
    print("Error: output directory " + args.output_dir + " already exists.")
    exit(-1)
    
with open(args.sxs_catalog_metadata, 'r') as file:
        sxs_catalog_metadata = json.load(file)
        file.close()
        
bbh_simulations = bbh_keys_from_simulation_keys(sxs_catalog_metadata['simulations'])

zen.download.matching("common-metadata.txt", "metadata.txt",
                      "rhOverM_Asymptotic_GeometricUnits_CoM.h5",
                      "Horizons.h5",
                      sxs_ids = bbh_simulations,
                      dry_run = args.dry_run,
                      highest_lev_only = args.highest_lev_only)

