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
  * `sxs_metadata_example.ipynb` -- how to download and manipulate zenodo metadata
  * `sxs_arxiv_metadata_example.ipynb` -- how to download and manipulate sxs metadata included in arXiv:1904.04831
  * `waveform_tutorial.ipynb` -- how to work with finite-radius waveforms and extrapolated waveforms, and how to compute the wave polarizations at a chosen sky location

Please send questions to questions@black-holes.org.

