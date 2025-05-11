import pandas as pd

# make MAG and metadata table
df_meta=pd.read_csv('skin_project_meta.csv')
df_meta_simple=df_meta.copy()
df_meta_simple['Site-Symmetry']=df_meta_simple['Site-Symmetry'].str.replace('-R','').str.replace('-L','').str.replace('-C','').str.replace('-null','').str.replace(':WGA','')
df_meta_simple['Disease']=df_meta_simple['Conditions'].apply(lambda x: x if x =='healthy' else 'disease')
df_meta_simple.loc[df_meta_simple['Site-Symmetry'].str.contains('Tn'), 'Site characteristic'] = 'foot'
df_meta_simple.loc[df_meta_simple['Site-Symmetry'].str.contains('Tw'), 'Site characteristic'] = 'foot'
df_meta_simple.loc[df_meta_simple['Site-Symmetry'].str.contains('Ph'), 'Site characteristic'] = 'foot'
df_meta_simple['Disease_detail']=df_meta_simple['Conditions'].apply(lambda x: x if x =='healthy' else 'disease')+'_'+df_meta_simple['Site characteristic']
df_meta_simple['Designation']=df_meta_simple['Country']+'_'+df_meta_simple['Designation'].astype("string")
df_meta_simple['Characteristic']=df_meta_simple['Site characteristic'].astype("string")
df_meta_simple=df_meta_simple.drop(['Age at time of sampling','Note','Date Collected','Total final reads','Method'],axis=1)
df_meta_simple=df_meta_simple.rename(columns={'Site characteristic':'Site_characteristic','Site-Symmetry':'Site_Symmetry','SAMPLE':'Sample'})
df_meta_simple=df_meta_simple.set_index('ID')
df_meta_simple=df_meta_simple.assign(Designation_gwasID=lambda x:x.Designation.str.replace('USA_','HV').str.replace('China_','China')).set_index('Designation_gwasID').drop_duplicates('Designation')[['Gender']]
samples=pd.read_csv('filtered_samples.txt',index_col=0,sep='\t',header=None).index.to_list()
df_meta_simple=df_meta_simple.loc[samples,:]
replace_to_numeric={'moist': 1,'sebaceous': 2,'dry': 3,'toenail':4,'scalp':5, 'healthy':0, 'dandruff':1, 'AD':2, 'USA':0, 'India':1, 'Singapore':2, 'China':3, 'M':1, 'F':2}
df_meta_simple=df_meta_simple.replace({"Gender":replace_to_numeric})
df_meta_simple=df_meta_simple.assign(FID=0).assign(IID=lambda x:x.index)[['FID','IID','Gender']]
df_meta_simple.columns=['#FID','IID','SEX']
df_meta_simple.to_csv('sex_info.txt',index=None,sep='\t')
