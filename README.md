# Supporting Data and Code for Flu Warning System

This repository contains the software and data associated with the manuscript "Multi-scale early warning system for influenza A spillovers" by Alfonsi T., Bernasconi A., Chiara M., and Ceri S.

The software takes as input files the metadata and HA genomes of Influenza viruses formatted as .xls files for metadata and .fasta files for the genomes. 

The software is organized in three stages, where each stage depends on the output of the previous one:
1. Sequence Alignment 
2. Data Processing
3. Data Visualization

Each stage is distributed as a separate docker container, which can be built and run with simple commands.

The following sections describe how to use the data available on GISAID with the software to reproduce the analysis described in the associated manuscript.

    > Disclaimers: 
    > - the data availability may change during time, therefore a few difference with the results described in the manuscript are expected.
    > - the Data Visualization stage includes pre-defined filters based on the Collection_Date of the input sequences. In order to fully display the output of the Data Processing stage, ensure that the downloaded sequences are collected between July 28, 2008 and January 10, 2010 for H1N1 data, and between January 1, 2019 and May 5, 2025 for H5N1 data. 

# Data download from GISAID
Download:
- isolates' metadata in gisaid_epiflu_isolates.xls
- isolates' HA genomes in gisaid_epiflu_sequence.fasta 

    > The fasta eader must have the format "DNA Accession no. | Isolate ID | Isolate name | Type | Lineage | Clade | Segment"
    > The date format must be 'YYYY-MM-DD'

# Software requirements

Download, install and run Docker Desktop on your computer. 

# Software setup

Using a terminal, enter each of the software directories and run the command `docker-compose build` to prepare the software stages.

# Software Usage

### 1. Sequence Alignment

Preprare the input for the `data_alignment` package by copying the isolates' HA genomes inside the proper directory:
- for H1N1 data: put the files into `data_alignment > inputs > h1n1`
- for H5N1 data: put the files into `data_alignment > inputs > h5n1`
Reference sequences for H1N1 and H5N1 are already provided in `data_alignment > inputs > references` and correspond to the following isolates/sequences on NCBI Nucleotide DB: 
- NC_026433.1 Influenza A virus (A/California/07/2009(H1N1)) segment 4 hemagglutinin (HA) gene, complete cds
- AF144305.1 Influenza A virus (A/Goose/Guangdong/1/96(H5N1)) hemagglutinin (HA) gene, complete cds

Open a terminal window, move inside the `data_alignment` folder, then separately launch the alignment procedure for H1N1 and H5N1 data with the commands:

- `docker-compose run data_alignment h1n1.sh` for aligning H1N1 files
- `docker-compose run data_alignment h5n1.sh` for aligning H5N1 files

At each command execution, the package will generate a pair of .fasta and .insertions.csv files in the `data_alignment > alignments` directory. 

### 2. Data Processing

The following instruction describes how to process the data concerining one serotype (H1N1 or H5N1). 

Prepare the input for the `data_processing` package by:
- copying the isolates' metadata and HA genomes in the directory `data_processing > input`. The software expects ONE .fasta file for the genomes and ONE .xls file for the metadata.
- copying the otuput of procedure *1. Sequence Alignment* in the directory `data_processing > alignments`. The software expects ONE .fasta file and ONE .insertions.csv file. 
- copy the relevant reference sequence from `data_alignment > inputs > references` into `data_processing > references > refseq.fasta`. The software expects ONE refseq.fasta file.

Open a terminal window, move inside the `data_processing` folder, then launch the software package with the command:
- `docker-compose run data_processing 1 1701` for H1N1
or 
- `docker-compose run data_processing 22 1728` for H5N1

    > The arguments (1 1701) or (22 1728) tell the software the HA coordinates range in the aligned sequences. 

At each command execution, the package will generate a `datawarehouse.sqlite` file. To avoid overwriting it on the next run, please rename this file as `h1n1.sqlite` or `h5n1.sqlite` according to the input. 

### 3. Data Visualization

Please copy the output of *2. Data Processing* (two .sqlite files contained in the `data_processing > output` directory) to the directory `data_visualization > data`. Ensure the files are named as:
- `h1n1.sqlite` for the output of *2. Warnings analysis* about H1N1
- `h5n1.sqlite` for the output of *2. Warnings analysis* about H5N1

Open a terminal window, move inside the `data_visualization` folder, then:
- launch the visualization server with the command `docker-compose up`
- open a browser window at the url (http://localhost:60119/)[http://localhost:60119/]

When necessary, shutdown the visualization server by pressing [Ctrl+C] on Windows / Linux and [Cmd+C] on MacOS. 

# License

See the LICENSE file. 