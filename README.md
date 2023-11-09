# Shareable-code

This repository will allow anyone in the GeoMIP community to share useful code to analyze data.
If you publish a paper based on GeoMIP data, you're welcome to publish your code here as well (or create a new repository to be turned into a DOI with Zenodo).

Feel free to use the discussion to request advices as well!

What's in the repo:

**read_data_model.m** A Matlab file used to analyse and produce some of the figures in https://acp.copernicus.org/articles/21/10039/2021/ (Identifying the sources of uncertainty in climate model simulations of solar radiation modification with the G6sulfur and G6solar Geoengineering Model Intercomparison Project (GeoMIP) simulations by Visioni et al., 2021). This allows, once the relevant CMIP6 data is downloaded, to read them all at once and produce some global mean figures.

**plot_maps.m** A Matlab file used to analyse and produce some of the figures in https://acp.copernicus.org/articles/21/10039/2021/ (Identifying the sources of uncertainty in climate model simulations of solar radiation modification with the G6sulfur and G6solar Geoengineering Model Intercomparison Project (GeoMIP) simulations by Visioni et al., 2021). Once read_data_model.m is used, this allows to make 2D map per each model and for the multi-model average, with strippling for model agreement/disagreement. Here used for temperature.
