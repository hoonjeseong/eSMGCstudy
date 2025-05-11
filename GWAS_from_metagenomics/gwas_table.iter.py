import sys
import pandas as pd
import numpy as np
import pingouin as pg
import pickle
from glob import glob
from Bio.SeqIO.FastaIO import SimpleFastaParser as sfp
from collections import defaultdict
from collections import Counter

# check human SNP PCA
dfPC=pd.read_csv('pca_result.eigenvec',index_col=1,sep='\t').drop('#FID',axis=1) #### <- add pca eignvector from PLINK2
dfPC=dfPC.assign(Group=lambda x:x.index.str[:2])
# make MAG and metadata table
df_meta=pd.read_csv('skin_project_meta.csv') #### <- add metadata
df_meta_simple=df_meta.copy()
df_meta_simple['Site-Symmetry']=df_meta_simple['Site-Symmetry'].str.replace('-R','').str.replace('-L','').str.replace('-C','').str.replace('-null','').str.replace(':WGA','')
df_meta_simple['Disease']=df_meta_simple['Conditions'].apply(lambda x: x if x =='healthy' else 'disease')
df_meta_simple.loc[df_meta_simple['Site-Symmetry'].str.contains('Tn'), 'Site characteristic'] = 'foot'
df_meta_simple.loc[df_meta_simple['Site-Symmetry'].str.contains('Tw'), 'Site characteristic'] = 'foot'
df_meta_simple.loc[df_meta_simple['Site-Symmetry'].str.contains('Ph'), 'Site characteristic'] = 'foot'
df_meta_simple['Disease_detail']=df_meta_simple['Conditions'].apply(lambda x: x if x =='healthy' else 'disease')+'_'+df_meta_simple['Site characteristic']
df_meta_simple['Designation_gwasID']=df_meta_simple['Country']+'_'+df_meta_simple['Designation'].astype("string")
df_meta_simple['Characteristic']=df_meta_simple['Site characteristic'].astype("string")
df_meta_simple=df_meta_simple.drop(['Age at time of sampling','Note','Date Collected','Total final reads','Method'],axis=1)
df_meta_simple=df_meta_simple.rename(columns={'Site characteristic':'Site_characteristic','Site-Symmetry':'Site_Symmetry','SAMPLE':'Sample'})

tmp=df_meta_simple.set_index('Sample')['ID'].to_dict()
df_MAG=pd.read_csv('MAG.abundance.txt',index_col=0,sep='\t') #### <- add MAG abundance table

df_meta_simple=df_meta_simple.set_index('ID')

# Function to perform CLR transformation
def clr_transformation(df):
    # Replace zeros with a small number to avoid division errors in log transformation
    df = df.replace(0, np.finfo(float).eps)
    # Compute the geometric mean of rows
    geometric_mean = df.apply(lambda x: np.exp(np.mean(np.log(x))), axis=1)
    # Perform CLR transformation
    clr_transformed = df.apply(lambda x: np.log(x / geometric_mean[x.name]), axis=1)
    return clr_transformed

# Applying the CLR transformation to the DataFrame
df_MAG_clr = clr_transformation(df_MAG.transpose()).transpose()

# filter MAG by mapped coverage more than 1X dp
# with mobilome(vMAG, pMAG) >= 30%, prokaryote(proMAG) >= 10% and fungi(eukMAG) >= 5% coverage\n",
df_breadth1=pd.read_csv('genome_breadth.1dp.txt',sep='\t').melt(id_vars='MAGs') #### <- add mapped coverage table
df_breadth1.columns=['MAGs','ID','value']
df_breadth_f1=(df_breadth1.loc[lambda x: ~x.ID.isnull()]
 .loc[lambda x:((x.MAGs.str.contains('_v_')) & (x.value>=0.3)) | ((x.MAGs.str.contains('_p_')) & (x.value>=0.3)) |\
      ((x.MAGs.str.contains('_euk_')) & (x.value>=0.05)) | ((x.MAGs.str.contains('_pro_')) & (x.value>=0.1))])

#select site for N Ch Gb Ac and pick random sample (with seed number)
df_for_test=df_meta_simple.loc[lambda x:((x.Country=='China')) | (x.Country=='Singapore') | (x.Country=='India') |((x.Site_Symmetry.isin(['Al','Ch','Gb','Ac'])) & (x.Country=='USA'))]
seed_number=int(sys.argv[1])
random_selected = df_for_test.groupby('Designation_gwasID').apply(lambda x: x.sample(1, random_state=seed_number)).drop('Designation_gwasID',axis=1).reset_index().set_index('ID')
random_selected=random_selected[['Designation_gwasID','Country','Gender','Site_characteristic','Conditions','Sample']]
df_MAG=df_MAG[random_selected.index.to_list()].T
df_MAG_clr=df_MAG_clr[random_selected.index.to_list()].T

#make prevalance table
ID_meta=(random_selected.set_index('Sample')['Designation_gwasID']
  .T.to_dict())
df_Prev=df_breadth_f1.loc[lambda x:x.newID.isin(ID_meta)].drop('ID',axis=1).pivot_table(columns='MAGs',values='value',index='newID')
df_Prev['ID']=df_Prev.index.map(ID_meta)
df_Prev=df_Prev.set_index('ID').fillna(0)
df_Prev[df_Prev > 0.0] = 1.0
# Filter columns where more than 50% of values are greater than 0.001
#threshold_value = 0.1; threshold_percentage = 0.2 #for the df_MAG
threshold_value = 0.0; threshold_percentage = 0.2 # for the df_Prev
df_Test=df_Prev.copy() #df_MAG or df_Prev
filtered_columns = [col for col in df_Test if (df_Test[col] > threshold_value).mean() > threshold_percentage]
df_Test = df_MAG_clr[filtered_columns]

# Greedy clustering based on R^2 > 0.9
correlation_matrix=df_Test.corr()
def greedy_clustering(correlation_matrix, threshold):
    clusters = []
    used = set()
    for column in correlation_matrix.columns:
        if column not in used:
            new_cluster = [column]
            used.add(column)
            # Find all columns that are not yet used and highly correlated with the current column
            for other_col in correlation_matrix.columns:
                if other_col not in used and correlation_matrix.at[column, other_col] ** 2 > threshold:
                    new_cluster.append(other_col)
                    used.add(other_col)
            clusters.append(new_cluster)
    return clusters
clusters = greedy_clustering(correlation_matrix, 0.9)
rep_set=[]
for i, cluster in enumerate(clusters):
    print(f"Cluster {i+1}: {cluster}")
    rep_set.append(cluster[-1])

##check confidency of pMAG and vMAG <- from JGI IMG site
IMGpr='IMGPR_plasmid_data.tsv'
IMGvr='IMGVR_all_Sequence_information-high_confidence.tsv'
prTaxon=defaultdict(set)
vrTaxon=defaultdict(set)
with open('pMAG_to_jgi.txt','r') as P: #blasting pMAG to JGI_PR
    next(P)
    for i in P:
        i=i.rstrip('\n').split('\t')
        t=[j.split('|')[0] for j in i[-1].split(',')]
        for j in t:
            prTaxon[i[0]].add(j)
with open('vMAG_to_jgi.txt','r') as P: #blasting vMAG to JGI_VR
    next(P)
    for i in P:
        i=i.rstrip('\n').split('\t')
        t=[j.split('|')[0] for j in i[-1].split(',')]
        for j in t:
            vrTaxon[i[0]].add(j)
for i in rep_set:
    if '_v_' in i:
        if not i in vrTaxon:
            rep_set.remove(i)
    elif '_p_' in i:
        if not i in prTaxon:
            rep_set.remove(i)

#order data frame
samples=pd.read_csv('filtered_samples.txt',index_col=0,sep='\t',header=None).index.to_list()
df_Test=df_Test[rep_set]
rename=random_selected.Designation_gwasID.to_dict()
df_Test=df_Test.assign(Sample=lambda x:x.index.map(rename)).set_index('Sample').loc[samples,:]
replace_to_numeric={'moist': 1,'sebaceous': 2,'dry': 3,'toenail':4,'scalp':5, 'healthy':1, 'dandruff':2, 'AD':3, 'USA':1, 'India':2, 'Singapore':3, 'China':4, 'M':1, 'F':2}
random_selected=random_selected.set_index('Designation_gwasID').loc[samples,:]
random_selected=random_selected.replace({"Country":replace_to_numeric, "Gender":replace_to_numeric, "Site_characteristic":replace_to_numeric, "Conditions":replace_to_numeric})

#save phenotyp and covariate txt
random_selected=random_selected.assign(FID=0).assign(IID=lambda x:x.index)[['FID','IID','Country','Gender','Site_characteristic','Conditions']]
random_selected.columns=['#FID','IID','Country','Gender','Site_characteristic','Conditions']
dfPC=pd.read_csv('pca_result.eigenvec',sep='\t',index_col=1).drop('#FID',axis=1)
pd.concat([random_selected,dfPC],axis=1).to_csv('iter'+str(sys.argv[1])+'/covariate.txt',sep='\t',index=None)

tmp=df_Test.columns.to_list()
df_Test=df_Test.assign(FID=0).assign(IID=lambda x:x.index)
df_Test=df_Test[['FID','IID']+tmp]
df_Test.columns=[['#FID','IID']+tmp]
df_Test.to_csv('iter'+str(sys.argv[1])+'/phenotype.txt',sep='\t',index=None)
