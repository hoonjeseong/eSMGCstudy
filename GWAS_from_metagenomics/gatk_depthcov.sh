#!/bin/bash

ls *.bqsr.bam > bam_list.txt
[ ! -d ../gatk_depth ] && mkdir ../gatk_depth

flag=0
while read line; do
    path=$line
    var=$var' '$line
    IFS='/' read -a FILENAME <<< "$line"
    if [[ $flag -eq 10 ]];then
        gatk DepthOfCoverage -O ../gatk_depth/${FILENAME[1]%.bqsr.bam}.cov -L ../interval.bed -R ../hg38/GCF_000001405.40_GRCh38.p14_genomic.fna -I ${path}
	flag=0
    else
	((flag++))
	gatk DepthOfCoverage -O ../gatk_depth/${FILENAME[1]%.bqsr.bam}.cov -L ../interval.bed -R ../hg38/GCF_000001405.40_GRCh38.p14_genomic.fna -I ${path} &
    fi
done < bam_list.txt 
