import sys
from collections import defaultdict

hit_ctg={}
blastn=sys.argv[1]
out=sys.argv[2]
taxon='NCBI_DB/NCBI_NT_20210518/nt.taxon'

with open(blastn,'r') as B:
  for i in B:
    i=i.rstrip('\n').split('\t')
    if float(i[2])<80:
      continue
    S=min(int(i[6]),int(i[7]))
    E=max(int(i[6]),int(i[7]))
    if not i[0] in hit_ctg:
      hit_ctg[i[0]]=defaultdict(list)
    if not 'len' in hit_ctg[i[0]]:
      hit_ctg[i[0]]['len']=int(i[-2])
    if i[1] in hit_ctg[i[0]]:
      tmp=[]
      for s,e in hit_ctg[i[0]][i[1]]:
        if s<=S<=e and E>=e:
          e=E
        elif s<=E<=e and S<=s:
          s=S
        elif s<=S<=e and s<=E<=e:
          s=s;e=e
        elif S<=s<=E and S<=e<=E:
          s=S;e=E
        tmp.append((s,e))
      hit_ctg[i[0]][i[1]]=tmp
    else:
      hit_ctg[i[0]][i[1]].append((S,E))

filter_set=defaultdict(dict)
for i in hit_ctg:
  for j in hit_ctg[i]:
    if not j=='len':
      size=sum([k[1]-k[0]+1 for k in hit_ctg[i][j]])
      if size/float(hit_ctg[i]['len']) >= 0.4:
        filter_set[i][j]=size/float(hit_ctg[i]['len'])

total_ref=set([j for i in filter_set for j in filter_set[i]])
flag=len(total_ref)
taxon_dic={}
tmp=set()
with open(taxon,'r') as T:
  for i in T:
    i=i.rstrip('\n').split('\t')
    if i[0] in total_ref:
      tmp.add(i[0])
      taxon_dic[i[0]]=i[1]
      flag-=1
    if flag==0:
      break

with open(out,'w') as O:
  O.write('Query\tTarget\tTaxon\tAlign_perc\n')
  for i in filter_set:
    for j in filter_set[i]:
      if j=='len':
        continue
      if not j in taxon_dic:
        continue
      O.write(i+'\t'+j+'\t'+taxon_dic[j]+'\t'+str(filter_set[i][j])+'\n')
