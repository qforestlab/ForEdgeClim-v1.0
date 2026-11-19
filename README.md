# ForEdgeClim-v1.0

**ForEdgeClim v1.0** is an R package for simulating three-dimensional forest microclimates at high spatial resolution (1 m³). The model resolves key energy balance components, including shortwave and longwave radiation, sensible and latent heat fluxes, and ground heat fluxes.

The model was originally developed to simulate microclimate gradients along **forest edge-to-core transects**, explicitly accounting for vertical and lateral radiative and thermal exchanges within forest canopies.

This repository contains the original **ForEdgeClim v1.0** model. A detailed description of the model structure, underlying physical equations, and evaluation is provided in the associated publication.

## Associated publication

Van de Walle et al. (2026). *ForEdgeClim v1.0: a 3D process-based microclimate model incorporating vertical and lateral radiative and thermal fluxes to simulate forest edge-to-core transitions*. Geoscientific Model Development.

https://gmd.copernicus.org/articles/19/4661/2026/

## Installation

Install the package directly from GitHub:

```r
install.packages("devtools")
devtools::install_github("EmmaVdW27/ForEdgeClim-v1.0")
```

## Getting started

An example demonstrating how to run **ForEdgeClim v1.0**, including the preparation of input data, model execution, and visualisation of simulation outputs, is provided in:

```text
Scripts/testing_foredgeclim.R
```

The example illustrates the simulation of forest microclimate conditions using terrestrial laser scanning (TLS) data and meteorological observations.
