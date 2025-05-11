#### Mapping
#read mapping to human genome
for i in fastq/*_1.fq.gz
do
  sample=${i#fastq/}
  sample=${sample%_1.fq.gz}
  bwa mem -t 16 hg38/GCF_000001405.40_GRCh38.p14_genomic.fna fastq/${sample}_1.fq.gz fastq/${sample}_2.fq.gz | samtools view -Shu -F4 - > bam/${sample}.bam
  samtools sort -O BAM -o bam/${sample}.sort.bam bam/${sample}.bam
done

#merge bamfiles by individual
mkdir merged_bam
cp make_merge.py bam
cd bam
python make_merge.py > merge_bam.sh
bash merge_bam.sh
cd ../
########################################

#### GATK
mkdir GATK && cd GATK
#1.CreateSequenceDictionary
gatk CreateSequenceDictionary -R ../hg38/GCF_000001405.40_GRCh38.p14_genomic.fna
#2.run GATK
for i in ../merged_bam/*bam
do
  Sample=${i#../merged_bam/}
  Sample=${i%.bam}
  gatk MarkDuplicates -I ${i} -M ${Sample}.dup.txt -O ${Sample}.sort_markdup.bam --READ_NAME_REGEX null
  java -jar ~/gatk/picard/build/libs/picard.jar AddOrReplaceReadGroups \
       I=${Sample}.sort_markdup.bam \
       O=${Sample}.markdup.RG.bam \
       RGID=${Sample} \
       RGLB=${Sample} \
       RGPL=ILLUMINA \ # change this option to the BGI for only China samples
       RGPU=${Sample} \
       RGSM=${Sample}
  gatk BaseRecalibrator -I ${Sample}.markdup.RG.bam -R ../hg38/GCF_000001405.40_GRCh38.p14_genomic.fna --known-sites ../hg38/GCF_000001405.40.gz -O ${Sample}.markdup_bqsr
  gatk ApplyBQSR -bqsr ${Sample}.markdup_bqsr -I ${Sample}.markdup.RG.bam -O ${Sample}.bqsr.bam
  gatk HaplotypeCaller -I ${Sample}.bqsr.bam -R ../hg38/GCF_000001405.40_GRCh38.p14_genomic.fna -O ${Sample}.vcf.gz -ERC GVCF
done
#3. calculate coverage
bash gatk_depthcov.sh
cd ../
########################################

#### GLIMPSE2
#1. glimpse imputation
ls ../GATK/*bqsr.bam > glimpse_input.txt
[ ! -d GLIMPSE_impute ] && mkdir -p GLIMPSE_impute
cd GLIMPSE_impute

flag=0
for i in {1..22} X
do
  REF=/DB/glimpse_db/1000G/db_sites/split/1000GP.chr${i}
  BAM=glimpse_input.txt
  IFS='/' read -a FILENAME <<< "$b"
  while IFS="" read -r LINE || [ -n "$LINE" ];
  do
    printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
    IRG=$(echo $LINE | cut -d" " -f3)
    ORG=$(echo $LINE | cut -d" " -f4)
    CHR=$(echo ${LINE} | cut -d" " -f2)
    REGS=$(echo ${IRG} | cut -d":" -f 2 | cut -d"-" -f1)
    REGE=$(echo ${IRG} | cut -d":" -f 2 | cut -d"-" -f2)
    OUT=imputed
    if [ $flag -eq 20 ];then
      GLIMPSE2_phase_static --threads 4 --bam-list ${BAM} --reference ${REF}_${CHR}_${REGS}_${REGE}.bin --output ${OUT}_${CHR}_${REGS}_${REGE}.bcf
      flag=0
      sleep 3m
    else
      flag=$(( flag + 1 ))
      GLIMPSE2_phase_static --threads 4 --bam-list ${BAM} --reference ${REF}_${CHR}_${REGS}_${REGE}.bin --output ${OUT}_${CHR}_${REGS}_${REGE}.bcf &
    fi
  done < /DB/glimpse_db/1000G/db_sites/chunks.chr${i}.txt
done

#2. filter samples by coverage > 0.1X
python filter_sample.py
bash run_bcf_filter.sh

#3. glimpse ligate
[ ! -d GLIMPSE_ligate ] && mkdir -p GLIMPSE_ligate
cd GLIMPSE_ligate
flag=0
for i in {1..22} X
do
  if [ $flag -eq 25 ];then
    LST=chr${i}.txt
    ls -1v ../GLIMPSE_impute/filter/imputed_chr${i}_*bcf > chr${i}.txt
    OUT=chr${i}_ligated.bcf
    GLIMPSE2_ligate_static --input ${LST} --output ${OUT} > chr${i}.log 2>&1
    flag=0
    sleep 3m
  else
    flag=$(( flag + 1 ))
    LST=chr${i}.txt
    ls -1v ../GLIMPSE_impute/filter/imputed_chr${i}_*bcf > chr${i}.txt
    OUT=chr${i}_ligated.bcf
    nohup GLIMPSE2_ligate_static --input ${LST} --output ${OUT} > chr${i}.log 2>&1 &
  fi
done

cat chr1.txt chr2.txt chr3.txt chr4.txt chr5.txt chr6.txt chr7.txt chr8.txt chr9.txt chr10.txt chr11.txt chr12.txt chr13.txt chr14.txt chr15.txt chr16.txt chr17.txt chr18.txt chr19.txt chr20.txt chr21.txt chr22.txt chrX.txt > total_chr.txt
GLIMPSE2_ligate_static --threads 72 --input total_chr.txt --output total_ligated.bcf
cd ..
########################################

#### run GWAS
[ ! -d GWAS ] && mkdir -p GWAS
cd GWAS
cp GLIMPSE_impute/total_ligated* ./
#1. filter bcf - MAF 5% and Max Missing 90%
bcftools view total_ligated.bcf --types snps |vcftools --vcf - --max-missing 0.9 --maf 0.05 --recode --stdout --recode-INFO-all |bgzip > total_ligated.vcf.gz
#2. pre-work for iterative GWAS #only for chromosome 1-22
python add_sex.py
plink2 --vcf total_ligated.vcf.gz --make-bed --out skin.human --update-sex sex_info.txt --chr 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22
plink2 --bfile skin.human --pca 10 --out pca_result
plink2 --bfile skin.human --make-rel --out distance_matrix
#3. iterative GWAS test
for i in $(seq 1 120);
do
  if [ -d iter$i ]; then
      continue
  fi
  mkdir iter${i}
  python gwas_table.iter.py $i > iter${i}/log
  plink2 --threads 24 --bfile ../skin.human --covar iter${i}/covariate.txt --covar-name Country,Gender,Conditions,PC1,PC2,PC3,PC4 --glm --pheno iter${i}/phenotype.txt --out iter${i}/gwas_results
done
#4. calculate mean P values across 100 tests
bash run_mean_value.total.sh
########################################
