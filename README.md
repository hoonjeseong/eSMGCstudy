This repository contains scripts and methods developed for a dual-purpose analysis of skin metagenomic sequencing data:

__- Human Genome Imputation via Skin Metagenomic samples__
A pipeline to extract and impute host (human) genomic information from skin metagenomic data. This method aligns human-origin reads to the reference genome, calls variants using GATK, imputes genotypes with GLIMPSE2, and supports downstream genome-wide association studies (GWAS) using PLINK2. The approach enables host genetic inference from non-invasive microbiome samples with substantial human DNA content.

__- BLASTN-Based Screening for False Positives in Mobile Contig Detection__
A post-processing method that applies BLASTN to assembled contigs to screen for false positives in viral & plasmid sequence identification. The method merges non-overlapping aligned regions per subject, applies a minimum query coverage threshold (≥40%). If viral matches are outvoted by non-viral domains, the contig is flagged as likely non-viral. This is designed to improve robustness in virus & plasmid detection from skin metagenomes.

_📄 the preprint: https://www.biorxiv.org/content/10.1101/2025.05.01.651599v1_

- Most scripts require path modifications to match your local environment, including references to programs and cluster-specific directories. Please contact us if you need further clarification or input file format descriptions: Dr. Hoon Je Seong (hoonje@ksnu.ac.kr)
