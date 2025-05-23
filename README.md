# Supporting Data and Code for Flu Warning System

This repository contains the software and data associated with the manuscript "Multi-scale early warning system for influenza A spillovers" by Alfonsi T., Bernasconi A., Chiara M., and Ceri S.

The software takes as input files the metadata and HA genomes of Influenza viruses formatted as .xls files for metadata and .fasta files for the genomes. 

The software is organized in three stages, where each stage depends on the output of the previous one:
1. Sequence Alignment 
2. Data Processing
3. Data Visualization

Each stage is distributed as a docker service to be executed independently.

The following sections describe how to use the data available on GISAID with the software to reproduce the analysis described in the associated manuscript.

# Data download from GISAID
Download:
- isolates' metadata in gisaid_epiflu_isolates.xls
- isolates' HA genomes in gisaid_epiflu_sequence.fasta 

    > Please note that: 
    > - The fasta header must have the format "DNA Accession no. | Isolate ID | Isolate name | Type | Lineage | Clade | Segment" and any date must be formatted as 'YYYY-MM-DD'.
    > - the data availability may change during time, therefore some difference with the results described in the manuscript are expected
    > - the *3. Data Visualization* stage includes a filter on the Collection_Date of the input sequences. In order to fully display your input, ensure the sequences are collected between July 28, 2008 and January 10, 2010 for H1N1 data, and between January 1, 2019 and May 5, 2025 for H5N1 data. 

### Reference sequences for aligning H1N1 and H5N1 genomes
Reference sequences for H1N1 and H5N1 are already provided in `inputs > references` and correspond to the following isolates/sequences on NCBI Nucleotide DB: 
- NC_026433.1 Influenza A virus (A/California/07/2009(H1N1)) segment 4 hemagglutinin (HA) gene, complete cds
- AF144305.1 Influenza A virus (A/Goose/Guangdong/1/96(H5N1)) hemagglutinin (HA) gene, complete cds

# Software requirements

Download, install and run Docker Desktop on your computer. 

# Software setup

Using a terminal, run in this folder the command `docker-compose build` to prepare the software stages.

# Software Usage

### 1. Sequence Alignment

Preprare the input for the `data_alignment` package by copying the isolates' HA genomes inside the proper directory:
- for H1N1 data: put the files into `inputs > H1N1`
- for H5N1 data: put the files into `inputs > H5N1`

The files should be named as gisaid_epiflu_isolates.xls and gisaid_epiflu_sequence.fasta.

Open a terminal window, then separately launch the alignment procedure for H1N1 and H5N1 data with the commands:

- `docker-compose down --remove-orphans && docker-compose run data_alignment H1N1`
or
- `docker-compose down --remove-orphans && docker-compose run data_alignment H5N1`

This stage terminates with the generation of two pairs (one for serotype) of .fasta and .insertions.csv files in the `alignments` directory.

### 2. Data Processing

The following instruction describes how to process the data concerining one serotype (H1N1 or H5N1). 

Open a terminal window, then separately launch the data processing stage for H1N1 and H5N1 data with the commands:
- `docker-compose down --remove-orphans && docker-compose run data_processing H1N1`
or 
- `docker-compose down --remove-orphans && docker-compose run data_processing H5N1`

This stage terminates with the generation of the files H1N1.sqlite and H5N1.sqlite in the `output` directory.

### 3. Data Visualization

Open a terminal window and:
- launch the visualization server with the command `docker-compose up --remove-orphans data_visualization`
- open a browser window at the url (http://localhost:60119/)[http://localhost:60119/]

When necessary, shutdown the visualization server by pressing [Ctrl+C] on Windows / Linux and [Cmd+C] on MacOS. 

# License

See the LICENSE file. 