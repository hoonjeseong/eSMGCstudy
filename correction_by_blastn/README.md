![image](https://github.com/user-attachments/assets/be04198f-ad9a-4163-97dd-fbb2c392b8a2)

----
#### __command__
'''
blastn -db NCBI_DB/NCBI_NT_20210518/nt -query [input:fasta] -out [output:blastn result] -num_threads [threads] -evalue 1e-10 -perc_identity 80 -outfmt '6 std qlen slen'
python filter_blastn.py [input: blastn result] [output: blastn coverage result]
'''
