# Sugarcane-MgCh

This repository contains the code used for the sequence analysis in this paper:

Eid, A., Mohan, C., Sanchez, S., Wang, D. and Altpeter, F.<BR>**Multiallelic, targeted mutagenesis of magnesium chelatase with CRISPR-Cas9 provides a rapidly scorable phenotype in highly polyploid sugarcane.**<BR>*Frontiers in Genome Editing*. 2021

## Usage

Basic usage is as follows:

```bash
$ scan.py [options] fastq...
```

Where options are:

Option | Description
-------|------------
  -o O | Use O as base name for output files (default: '{}').
  -p   | Output percentages where appropriate instead of counts.

Fastq files can be compressed with gzip. 

The program generates three sets of tables, each one composed of several sheets. The first sheet
is a description of the contents of the other sheets. Files are named according to the following pattern:

```
  O.tableN.sheetM.txt
```

where O is set with the -o option, N is the table number, and M is the sheet number. All files are 
tab-delimited.

In addition, the program also generates one matrix for each input file, containing the base frequencies for
all positions in sgRNA1. These files are named as follows:

```
  O.matN.txt
```

where N is the index of the input file.


