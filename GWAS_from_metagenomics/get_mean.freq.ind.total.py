from glob import glob
import pandas as pd
import numpy as np
from collections import defaultdict
import pickle
import sys
mID=sys.argv[1]
SNPsig=[]
#with open('/BiO/Live/s4645/Project/Skin/human/GLIMPSE_ligate/total/maf5/not_ld_prune/no_Nreads/pickle_total/'+mID, 'rb') as handle:
with open(mID, 'rb') as handle:
    ID=mID.split('/')[-1].split('.')[0]
    gID=ID.split('/')[-1].split('_')[1]
    print (ID)
    b = pickle.load(handle)
    if b:
        for j in b:
            if len(b[j]['P'])<90:
                print (ID,len(b[j]['P']))
                break
            SNPsig.append([ID, gID, j, len([i for i in b[j]['P'] if i<5e-8]), len(b[j]['P']), np.mean(b[j]['P']), np.median(b[j]['P'])])
SNPsig=pd.DataFrame(SNPsig)
SNPsig.columns=['ID','gID','Pos','nSig','nTest','meanP','medP']
SNPsig.to_csv('pickle_freq/'+ID+'.freq.csv',sep='\t')
