# catalog_tools
Example scripts to interact with the SXS Catalog (https://black-holes.org/waveforms) and its data, hosted on Zenodo (https://zenodo.org). The data in the catalog is described in https://arxiv.org/abs/1904.04831.

This repo contains the following python scripts:
  * `get_sxs_bbh_catalog.py` -- download waveforms, horizon quantities, and metadata for 
  all binary-black-holes simulations in the SXS catalog
  * `get_sxs_public_metadata.py` -- download publicly available metadata for the catalog into JSON files
  * `convert_sxs_to_lvc.py` -- convert an SXS waveform into the LVC format (BETA)
  * `convert_sxs_metdata_txt_to_json.py` -- convert an SXS metadata.txt file to JSON format
  * `compare_sxs_converted_lvc.py` -- compare two different SXS waveforms in the LVC format, checking if attributes, splines, and time series are the same
  * `compare_sxs_vs_lvc.py` -- compare SXS data in the LVC format to data in the SXS format, checking that they agree

This repo contains the following notebooks in `Examples/`:
  * `sxs_catalog_download_example.ipynb` -- how to download and interact with catalog data
  * `sxs_metadata_example.ipynb` -- how to download and manipulate sxs metadata provided at data.black-holes.org
  * `sxs_arxiv_metadata_example.ipynb` -- how to download and manipulate sxs metadata included in arXiv:1904.04831
  * `sxs_zenodo_metadata_example.ipynb` -- how to download and manipulate zenodo metadata
  * `waveform_tutorial.ipynb` -- how to work with finite-radius waveforms and extrapolated waveforms, and how to compute the wave polarizations at a chosen sky location

## Dependencies
  * All scripts require python 3 (https://python.org) and have been tested using python 3.7.
  * `get_sxs_bbh_catalog.py` and `get_sxs_public_metadata.py` depnd on the `sxs` python module (https://github.com/sxs-collaboration/sxs), which requires python 3.
  * `convert_sxs_to_lvc.py`, `compare_sxs_converted_lvc.py`, `compare_sxs_converted_lvc.py`, and `compare_sxs_vs_lvc.py` depend on `romspline` (https://github.com/crgalley/romspline), which is now compatible with python 3.
  * To install the prerequisites, you can run the command `pip install sxs romspline`, or add the text `!pip install sxs romspline` to the beginning of a jupyter notebook.

## Questions

Please send questions to questions@black-holes.org.

