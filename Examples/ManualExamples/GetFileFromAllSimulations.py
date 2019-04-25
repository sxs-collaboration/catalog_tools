# This script makes a directory for each resolution of each simulation 
# in the SXS-collaboration catalog

import json
import argparse
import os
import requests

parser = argparse.ArgumentParser(description="Makes a directory hierarchy for each resolution of each simulation in the SXS collaboration catalog")

parser.add_argument('--res_input', help='Path to JSON file sxs_catalog_public_resolutions_available.json')
parser.add_argument('--zenodo_input', help='Path to JSON file sxs_catalog_public.json')
parser.add_argument('--catalog_root', help='Path to root directory of folder hierarchy for catalog data (via SetupCatalogDirectories.py)')
parser.add_argument('--file', help='name of file to download')

args = parser.parse_args()

# Read in a JSON file whose keys are the SXS ID numbers (SXS:BBH:####)
# for simulations in the catalog, and whose values are lists of 
# resolutions available (named as integers).

with open(args.res_input, 'r') as file:
    resolutions = json.load(file)

with open(args.zenodo_input, 'r') as file:
    catalog_json = json.load(file)

# Loop over simulations.
root_dir_name = args.catalog_root

for simulation in resolutions:
    print(simulation)

    simulation_dir_name = simulation.replace(':', '_')
    simulation_path = root_dir_name + "/" + simulation_dir_name
    simulation_json = catalog_json[simulation]

    simulation_files = simulation_json["files"]
    
    for res in resolutions[simulation]:
        for file in simulation_files:
            filename_split = file['filename'].split('/')
            if (str(filename_split[-1]) == args.file):
                if (str(filename_split[-2].split('Lev')[-1]) == str(res)):
                    url = file['links']['download']
                    output_filename = simulation_path + "/Res" + str(res) + "/" + args.file
                    print(url + " -> " + output_filename)
                    download = requests.get(url, allow_redirects=True)
                    open(output_filename, 'wb').write(download.content)
