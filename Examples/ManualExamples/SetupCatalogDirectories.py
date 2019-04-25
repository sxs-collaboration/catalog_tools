# This script makes a directory for each resolution of each simulation 
# in the SXS-collaboration catalog

import json
import argparse
import os

parser = argparse.ArgumentParser(description="Makes a directory hierarchy for each resolution of each simulation in the SXS collaboration catalog")

parser.add_argument('--input', help='Path to JSON file sxs_catalog_public_resolutions_available.json')
parser.add_argument('--output', help='Path to a directory where the directory hierarchy will be created')

args = parser.parse_args()

# Read in a JSON file whose keys are the SXS ID numbers (SXS:BBH:####)
# for simulations in the catalog, and whose values are lists of 
# resolutions available (named as integers).

with open(args.input, 'r') as file:
    resolutions = json.load(file)

# Loop over simulations.
root_dir_name = args.output
os.mkdir(root_dir_name)
for simulation in resolutions:
    simulation_dir_name = simulation.replace(':', '_')
    simulation_path = root_dir_name + "/" + simulation_dir_name
    os.mkdir(simulation_path)
    for res in resolutions[simulation]:
      os.mkdir(simulation_path + "/Res" + str(res))
