#### __Human Genome Imputation via Metagenomics__

This repository includes a method for human genome imputation from skin metagenomic sequencing data, which typically contains a substantial proportion of host (human) reads. 
The pipeline extracts these human-aligned reads, aligns them to the reference genome, and performs genotype imputation and association analysis.

__Key steps:__

_1. Read Mapping: Human reads are aligned to the GRCh38 reference genome using BWA-MEM._

_2. Variant Calling: The GATK Best Practices pipeline is followed, including duplicate marking, read group assignment, base quality recalibration, and variant calling via HaplotypeCaller in GVCF mode._

   _⚠️ For samples sequenced in China, the platform identifier RGPL should be set to BGI instead of ILLUMINA._

_3. Genotype Imputation: Low-coverage variant data is imputed using GLIMPSE2, suitable for metagenomic datasets._

_4. Iterative GWAS: Genome-wide association studies are conducted iteratively using PLINK2._

This approach enables host genome inference from microbiome data in settings where dedicated host sequencing is unavailable or not feasible.

📁 Please check the cmd.sh file for a step-by-step implementation.
📄 For detailed methodology and results, see the preprint: https://www.biorxiv.org/content/10.1101/2025.05.01.651599v1
