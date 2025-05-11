from glob import glob
import gzip as gz
from collections import defaultdict
import numpy as np
covD=defaultdict(list)
for c in glob('*_stats_coverage.txt.gz'):
  flag=0
  with gz.open(c,'r') as C:
      for i in C:
          i = i.decode('utf-8').split()
          flag+=1
          if flag<=2:
              continue
          if i[2]:
              covD[i[1]].append(float(i[2]))

filtered_list=[]
countryD=defaultdict(int)
for i in covD:
    if np.mean(covD[i])<0.1:
        continue
    filtered_list.append(i)
    if 'China' in i:
        countryD['CN']+=1
    elif 'HV' in i:
        countryD['US']+=1
    elif 'Singapore' in i:
        countryD['SG']+=1
    elif 'India' in i:
        countryD['ID']+=1
    elif 'Italy' in i:
        countryD['IT']+=1

with open('filtered_samples.txt','w') as F:
    for i in filtered_list:
        if 'China' in i:
            country='CN'
        elif 'HV' in i:
            country='US'
        elif 'Singapore' in i:
            country='SG'
        elif 'India' in i:
            country='ID'
        elif 'Italy' in i:
            country='IT'
        F.write(i+'\n')
