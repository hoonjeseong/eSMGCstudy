flag=0
for i in iter1/gwas*linear;
do
  i=${i#iter1/}
  ID=${i#gwas_results.}
  ID=${ID%.glm.linear}
  if [ -f pickle_freq/$ID.freq.csv ]; then
  #if [ -f pickle_total/$ID.pickle ]; then
      continue
  fi
  echo $ID
  if [ $flag == 50 ]; then
      python get_mean.ind.total.py $i > pickle_total/$i.log
      python get_mean.freq.ind.total.py pickle_total/$ID.pickle > pickle_freq/$ID.log
      sleep 3m
      flag=0
  else
      python get_mean.ind.total.py $i > pickle_total/$i.log &
      python get_mean.freq.ind.total.py pickle_total/$ID.pickle > pickle_freq/$ID.log &
      flag=$((flag+1))
  fi
done
