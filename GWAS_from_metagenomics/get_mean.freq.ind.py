from glob import glob
from collections import defaultdict
import pickle
import sys
mID=sys.argv[1]
SNPs=defaultdict(int)

for i in glob('iter*/'+mID):
    ID=i.split('/')[0]
    try:
        seed_number=int(i.split('/')[0].lstrip('iter'))
        if seed_number<121:
            #CHROM	POS	ID	REF	ALT	PROVISIONAL_REF?	A1	OMITTED	A1_FREQ	TEST	OBS_CT	BETA	SE	T_STAT	P	ERRCODE
            #1	205546793	1:205546793:C:T	C	T	Y	C	T	0.431193	ADD	327	-1.4666	0.255448	-5.74128	2.19963e-08	.
            with open(i,'r') as S:
                next(S)
                for l in S:
                    l=l.rstrip('\n').split('\t')
                    if l[9]!='ADD' or l[15]!='.' or float(l[14])> 5e-8: #3.68e-10:
                        continue
                    SNPs[l[2]]+=1
    except:
        continue
with open('pickle_freq/'+mID.split('.')[1]+'.freq.pickle', 'wb') as handle:
    pickle.dump(SNPs, handle, protocol=pickle.HIGHEST_PROTOCOL)
