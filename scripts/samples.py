import os
import pandas as pd
import re

aneuploidy_home = "/scratch/chd5n/aneuploidy/"
aneuploidy_raw_annotations = aneuploidy_home + "raw-data/annotations/"

# read all samples
os.chdir(aneuploidy_raw_annotations)
dat = pd.read_csv("gdc_sample_sheet.2019-07-13.tsv", sep="\t")

# explore all samples
dat.shape
dat.columns
dat['Case ID'].describe()
dat['Sample Type'].value_counts()

# select normal samples (blood and adjacent tissue)
normal_tissue = ['Blood Derived Normal', 'Solid Tissue Normal']
r = dat['Sample Type'].isin(normal_tissue)
dat_norm = dat.loc[r,:]
dat_norm = dat_norm.reset_index(drop=True)

# explore normal samples
dat_norm.shape
dat_norm['Case ID'].describe()
dat_norm['Sample Type'].value_counts()

# select duplicate normal samples
instances = dat_norm['Case ID'].value_counts()
duplicates = instances[instances>1]
dat_norm_dup = pd.concat([dat_norm.loc[dat_norm['Case ID']==i,:] for i in duplicates.index], ignore_index=True)

# explore duplicate normal samples
dat_norm_dup.shape
dat_norm_dup['Case ID'].describe()
dat_norm_dup['Sample Type'].value_counts()
instances.value_counts()

# write duplicate normal samples
dat_norm_dup.to_csv("samples-normal-duplicate_2019-07-13.tsv", sep="\t")

# select unique normal samples (exclude duplicates for now)
uniques = instances[instances==1]
dat_norm_unq = pd.concat([dat_norm.loc[dat_norm['Case ID']==i,:] for i in uniques.index], ignore_index=True)

# explore unique normal samples
dat_norm_unq.shape
dat_norm_unq['Case ID'].describe()
dat_norm_unq['Sample Type'].value_counts()

# write unique normal samples
dat_norm_unq.to_csv("samples-normal-unique_2019-07-13.tsv", sep="\t")

# select tumor samples
r = dat['Sample Type'] == 'Primary Tumor'
dat_tum = dat.loc[r,:]
dat_tum = dat_tum.reset_index(drop=True)

# explore tumor samples
dat_tum.shape
dat_tum['Case ID'].describe()
dat_tum['Sample Type'].value_counts()

# select duplicate tumor samples
instances = dat_tum['Case ID'].value_counts()
duplicates = instances[instances>1]
dat_tum_dup = pd.concat([dat_tum.loc[dat_tum['Case ID']==i,:] for i in duplicates.index], ignore_index=True)

# explore duplicate tumor samples
dat_tum_dup.shape
dat_tum_dup['Case ID'].describe()
instances.value_counts()

# write duplicate tumor samples
dat_tum_dup.to_csv("samples-tumor-duplicate_2019-07-13.tsv", sep="\t")

# select unique tumor samples (exclude duplicates for now)
uniques = instances[instances==1]
dat_tum_unq = pd.concat([dat_tum.loc[dat_tum['Case ID']==i,:] for i in uniques.index], ignore_index=True)

# explore unique tumor samples
dat_tum_unq.shape
dat_tum_unq['Case ID'].describe()
dat_tum_unq['Sample Type'].value_counts()

# write unique tumor samples
dat_tum_unq.to_csv("samples-tumor-unique_2019-07-13.tsv", sep="\t")

# build a pilot set of five subjects from the intersection of unique tumor and normal samples
    # dat_norm_dup
    # dat_norm_unq
    # dat_tum_dup
    # dat_tum_unq
unq_set = pd.merge(dat_norm_unq, dat_tum_unq, how='inner', on=['Case ID', 'Data Category', 'Data Type', 'Project ID'], suffixes=('_norm','_tum'))
header = ['Case ID', 'Data Category', 'Data Type', 'Project ID', 'File ID_norm', 'File Name_norm', 'Sample ID_norm', 'Sample Type_norm', 'File ID_tum', 'File Name_tum', 'Sample ID_tum', 'Sample Type_tum']
unq_set = unq_set[header]
pilot_set = unq_set.iloc[0:5,:]
# pilot set for slurm-based download (has never worked before)
dn_pilot_set = unq_set.iloc[5:10,:]

# write pilot set
pilot_set.to_csv("pilot-set_2019-07-13.tsv", sep="\t", index=False)

# build pilot manifest

## select file id and file name for norm and tum
norm_fields = ['File ID_norm', 'File Name_norm']
tum_fields = ['File ID_tum', 'File Name_tum']

df_norm = pilot_set[norm_fields]
field_dict = {'File ID_norm': 'id', 'File Name_norm': 'filename'}
df_norm = df_norm.rename(columns=field_dict)

df_tum = pilot_set[tum_fields]
field_dict = {'File ID_tum': 'id', 'File Name_tum': 'filename'}
df_tum = df_tum.rename(columns=field_dict)

## concatenate file ids and file names
pilot_id_name = pd.concat([df_norm, df_tum], ignore_index=True)

## intersect with manifest
full_manifest = pd.read_csv("gdc_manifest_20190713_183544.txt", sep="\t")
pilot_manifest = pd.merge(pilot_id_name, full_manifest, how='inner', on=['id', 'filename'])

# write pilot manifest
pilot_manifest.to_csv("pilot-manifest_2019-07-13.tsv", sep="\t", index=False)
# pilot set for slurm-based download
pilot_manifest.to_csv("dn-pilot-manifest_2019-07-13.tsv", sep="\t", index=False)



# explore hg19 samples
keep_rows = []
for i in dat_norm.index:
    if (len(re.findall('hg19',dat_norm.loc[i,'File Name']))>0):
        keep_rows.append(i)

dat_norm_hg19 = dat_norm.iloc[keep_rows,:]
dat_norm_hg19 = dat_norm_hg19.reset_index(drop=True)
