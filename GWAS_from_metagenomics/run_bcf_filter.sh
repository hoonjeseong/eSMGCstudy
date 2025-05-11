mkdir filter
for i in *.bcf
do
  bcftools view -S filtered_samples.txt $i -o filter/$i -O u
  bcftools index -c filter/$i
done
