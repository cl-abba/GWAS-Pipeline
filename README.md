# GWAS-pipeline

This repository describes steps in genome-wide association studies (GWAS) and post-GWAS analyses.  Compliment to the methods paper, "From genetic association to forensic prediction: Computational methods and tools for identifying phenotypically informative SNPs" which represents the second chapter of my PhD thesis.

## Introduction

This example uses vcf files from over 5,000 individuals from the Canadian Partnership for Tomorrow's Health (CanPath) dataset and investigates genetic associations with categorical eye colours (blue, green, hazel and brown). Categorical eye colour phenotype was self-reported by participants, along with province, age and sex. The command line scripts (.sh) are formatted for running on a server with a scheduler, such as SLURM. Some may require program installation if the server is not already equipped with the necessary programs. URLs are provided for each program utilized. It is reccomended to read the documentation before using a program for the first time, this will provide an understanding of the flags or parameters, how to set up the proper input files, what output files will look like and how to interpret them. Figures for the R visualization scripts are provided.

## Materials and methods

To replicate this protocol the following materials are required:

1. Phenotype file
2. Imputed genotype file (or files if separated by chromosome)
3. This example utilizes both imputed and not imputed genotype files, but the protocol still applies if only imputed files are available
4. Storage and workspace that can be used to storing output files, running jobs and downloading programs

The methodology is as follows:

1. Data storage, handling and cleaning (if necessary) - PLINK, bcftools, R
2. Pruning and quality control (QC) - PLINK, bcftools, vcftools
3. GWAS - PLINK, SAIGE
4. Meta-analysis - Metasoft
5. Visualization - R
6. Postliminary analyses (COJO and FINEMAP examples) - GCTA, LDSTORE2, FINEMAP
