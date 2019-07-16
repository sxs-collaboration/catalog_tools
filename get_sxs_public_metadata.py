# Distributed under the MIT License.
# See LICENSE.txt for details.

import argparse
import datetime
import numpy as np
import json
import os
import requests
import sxs
from sxs import zenodo as zen
import sys


def bbh_keys_from_simulation_keys(simulation_keys):
    bbh_simulations = []
    for simulation in simulation_keys:
        if simulation.split(':')[-2] == "BBH":
            bbh_simulations.append(simulation)
    return bbh_simulations


def resolutions_for_simulation(sxs_id, sxs_catalog):
    """Returns a list of the available resolutions for a given sxs_id
    and sxs catalog metadata"""
    resolutions = []
    url = sxs_catalog['simulations'][sxs_id]['url']
    files = sxs_catalog['records'][url]['files']
    for file in files:
        split_filename = file['filename'].split('/')
        if (str(split_filename[-1]) == "Horizons.h5"):
            resolutions.append(int(split_filename[-2].split('Lev')[-1]))
    return sorted(resolutions)


p = argparse.ArgumentParser(description="Get SXS public BBH metadata")
p.add_argument("--output_dir",
               help="Name of output directory (must not exist)",
               required=True)
args = p.parse_args()

try:
    os.mkdir(args.output_dir)
except:
    print("Error: output directory " + args.output_dir + " already exists.")
    exit(-1)

# Get current black-holes.org catalog.json from black-holes.org
print("Getting " + args.output_dir + "/sxs_catalog.json")
request = requests.get("https://data.black-holes.org/catalog.json",
                       headers={'accept': 'application/citeproc+json'})
sxs_catalog = request.json()

# Save the metadata
with open(args.output_dir + "/sxs_catalog.json", 'w') as file:
    file.write(json.dumps(sxs_catalog))

# Generate sxs_catalog_resolutions.json
simulation_keys = sxs_catalog['simulations'].keys()
bbh_simulations = bbh_keys_from_simulation_keys(simulation_keys)
print("Generating " + args.output_dir + "/sxs_catalog_resolutions.json")
resolutions_available = {}
for simulation in bbh_simulations:
    resolutions_available[simulation] = \
        resolutions_for_simulation(simulation, sxs_catalog)

with open(args.output_dir + "/sxs_catalog_resolutions.json", 'w') as file:
    file.write(json.dumps(resolutions_available))
