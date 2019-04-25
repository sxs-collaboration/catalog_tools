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


def resolutions_for_simulation(sxs_id, zenodo_metadata):
    """Returns a list of the available resolutions for a given sxs_id
    and zenodo metadata metadata"""
    resolutions = []
    files = zenodo_metadata[sxs_id]['files']
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

# Get zenodo metadata: sxs_zenodo.json
print("SXS zenodo metadata -> " + args.output_dir + "/sxs_zenodo.json")
md = zen.api.Records.search(q='communities:sxs AND access_right:open')
catalog_json = {}
for simulation in md:
    sxs_id = simulation['title'].split(' ')[-1]
    if (len(sxs_id.split(':')) == 3):
        if (sxs_id.split(':')[-2] == 'BBH'):
            catalog_json[sxs_id] = simulation
with open(args.output_dir + "/sxs_zenodo.json", 'w') as file:
    file.write(json.dumps(catalog_json))

# Generate sxs_zenodo_resolutions.json
print("Generating " + args.output_dir + "/sxs_zenodo_resolutions.json")
resolutions_available = {}
for simulation in catalog_json:
    resolutions_available[simulation] = \
        resolutions_for_simulation(simulation, catalog_json)

with open(args.output_dir + "/sxs_zenodo_resolutions.json", 'w') as file:
    file.write(json.dumps(resolutions_available))

# Get arxiv.org catalog.json
print("Getting " + args.output_dir + "/sxs_catalog_arxiv.json")
request = requests.get("https://arxiv.org/src/1904.04831/anc/sxs_catalog.json",
                       headers={'accept': 'application/citeproc+json'})
sxs_catalog_arxiv = request.json()
with open(args.output_dir + "/sxs_zenodo_arxiv.json", 'w') as file:
    file.write(json.dumps(sxs_catalog_arxiv))
