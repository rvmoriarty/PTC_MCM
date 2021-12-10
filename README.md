# Therapeutic Vaccination
Scripts used for analyzing SIVmac239M barcode sequences in therapeutically vaccinated MCM. THIS IS A WORK IN PROGRESS. 

This data is generated using a series of R scripts to demultiplex FASTQ files and identify SIVmac239M barcodes present in the plasma longitudinally and in lymph nodes at necropsy. 

## Scripts used 
The first two scripts used to demultiplex data and identify barcodes are provided to us courtesy of the Keele Lab (NCI Frederick). We change one script to reflect where our data is, and this one script calls the other. This results in a multi-page excel file. 

### Barcode_analysis.R
This is a work-in-progress R script that reads the excel files, extracts data, calculates diversity indices, filters based on frequency, and returns a reformatted file containing the animal number, days post infection, replicate number, number of barcodes identified per sample in total, number of barcodes identified over 1%, 3%, and 5% of the population, frequency of each barcode, and diversity of the sample.

### Additional scripts are in progress to identify the relative frequencies of non-barcoded SIVmac239 and barcoded SIVmac239M in samples following re-challenge. 
