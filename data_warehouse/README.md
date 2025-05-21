This repository includes the website/visualization tool used for illustrating the results discussed in the manuscript "Multi-scale early spillover warnings for influenza A viruses" by Alfonsi T., Bernasconi A., Chiara M., and Ceri S.

This document is organized as follows:
- In the _Data Requirements_ section, we describe the data used by this software.
- The _Software Requirements_ section lists the packages needed for running the website.
- The section _Usage Guide_ explains how to launch the website locally.

# Data Requirements

The website expects the following files to be available:
- `data/h1n1_08-09.sqlite`
- `data/h5n1_2019-2025-01-27.sqlite`
- `data/sensitivity_specificity10.parquet`

in order to collect these files, refer to documentation provided about the Data Processing part of this work (see file `data_processing/README.md`).

# Software Requirements

This website is shipped within a Docker container, therefore all the software dependencies are handled by Docker itself. The only requirement to run the Docker container, is to download and launch the latest version of Docker Desktop. 

# Usage Guide

To start the website, open a terminal inside the `data_warehouse` directory and run the command `docker-compose build && docker-compose down && docker-compose up;`.