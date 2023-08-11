# Characterizing SIV populations in SIVmac239M-infected MCM that exhibited post-treatment control
Scripts used for analyzing SIVmac239M barcode sequences, the proportion of rechallenge virus present, and identifying SNPs and CTL epitope variants in MCM. 

This data is generated using a series of R scripts to demultiplex FASTQ files and identify SIVmac239M barcodes present in the plasma longitudinally. 

## Scripts used 
The first two scripts used to demultiplex data and identify barcodes are provided to us courtesy of the Keele Lab (NCI Frederick). We change one script to reflect where our data is, and this one script calls the other. This results in a multi-page excel file. 

### 1. SIVmac239M_barcode_analysis.R
This R script reads the excel files generated by the Keele Lab R scripts, extracts and filters barcode data, calculates barcode diversity, converts barcode frequencies to lineage-specific viral load, and performs computational simulations to estimate the proportion of rebounding lineages per pre-ART viral load bin and logistic regressions for estimating the impact of transient viral replication on the composition of viral lineages post-depletion. 

Inputs required: 
- TSV of known barcodes
- CSV of animal and sample viral loads
- Excel sheets generated by the Keele Lab barcode identification R script
- Manually filtered CSV of barcodes (I didn't want to spend the time to code this in)

Outputs generated: 
- CSV files containing the barcodes present over a given threshold for each animal in copies/mL (used to make the manually adjusted CSVs)
- CSV of rebound statistics calculated by the computational simulation of rebounding lineages by pre-ART VL

Additional outputs that are not written to a file can be found when running the script. 

### 2. WholeGenomeVariantAnalysis.ipynb 
This Jupyter notebook is designed to take whole-genome sequence data from paired-end FASTQ files, merge reads, map to a reference,  and identify epitope variants and SNPs. 

Python packages required: 
- samtools (version 1.9)
- os
- pathlib
- tempfile
- shutil
- subprocess
- quasitools
- pandas
- glob
- pyfastx
- re
- collections
- Bio.Seq
- matplotlib (pyplot)

Additional programs required: 
- SNPGenie and its dependencies
- SNPEff and its dependencies 
- BBtools
- Short amplicon normalization scripts (See ref) 
  
Inputs required:
- GTF and .fasta file for your reference
- CSV file of epitope sequences
- List of illumina indices 

Outputs generated: 
- TXT and CSV files for each animal containing epitope sequences
- Whole genome merged .fastq.gz files and .bam files, and consensus sequences for gag, env, and pol
- VCF file containing whole-genome variants (and one annotated by SNPEff) 
- SNPGenie folder containing data for each  sample

  ### WholeGenomeEpitopeAnalysis.R

  This R script is used to  organize the CSV/TXT files
 
