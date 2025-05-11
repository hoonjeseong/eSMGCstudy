![image](https://github.com/user-attachments/assets/be04198f-ad9a-4163-97dd-fbb2c392b8a2)

----
This method utilizes the output of BLASTN alignments against the NCBI NT database to classify assembled contigs. It addresses the common issue of multiple partial alignments within a single contig by evaluating each alignment's coverage and taxonomic assignment. For each subject, non-overlapping aligned regions on the query contig are merged to compute cumulative coverage.

----
#### __command__
1. run blastn with NCBI nt db
   
`blastn -db NCBI_DB/NCBI_NT_20210518/nt -query [input:fasta] -out [output:blastn result] -num_threads [threads] -evalue 1e-10 -perc_identity 80 -outfmt '6 std qlen slen'`

2. filter coverage over 40% with 80% identity
   
`python filter_blastn.py [input: blastn result] [output: blastn coverage result]`
