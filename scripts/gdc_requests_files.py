#!/usr/bin/env python

    # this script queries the gdc archive via the search and retrieve api and
    # returns the files_res object (results from files endpoint query)

import io
import json
import os
import pandas as pd
import requests
import re

data_home = '/scratch/chd5n/aneuploidy/'
anno_path = data_home + 'raw-data/annotations/'

# request WXS file_id entities from TCGA-COAD and READ
files_endpt = 'https://api.gdc.cancer.gov/files/'

fields = [
    'file_name',
    'file_id',
    'file_size',
    'md5sum',
    'state',
    'submitter_id',
    'experimental_strategy', # 'WXS'
    'data_category', # 'Sequencing Reads'
    'data_format', # 'BAM'
    'cases.case_id',
    'cases.project.project_id',
    'cases.primary_site',
    'cases.demographic.gender',
    'cases.demographic.race',
    'cases.demographic.year_of_birth',
    'cases.demographic.vital_status',
    'cases.demographic.days_to_death',
    'cases.demographic.age_at_index',
    'cases.samples.sample_type', # tumor vs normal
    'cases.samples.is_ffpe',
    'cases.samples.sample_id',
    'cases.tissue_source_site.name',
    'cases.exposures.weight',
    'cases.exposures.bmi',
    'cases.exposures.height'

    ## nested nature makes more suitable for special lookup
    #'analysis.metadata.read_groups.read_group_qcs.total_sequences',

    ## nested nature makes more suitable for invididual lookup
    # 'analysis.metadata.read_groups.target_capture_kit_name',
    # 'analysis.metadata.read_groups.target_capture_kit',
    # 'analysis.metadata.read_groups.target_capture_kit_target_region',
    # 'analysis.metadata.read_groups.target_capture_kit_catalog_number',
    # 'analysis.metadata.read_groups.target_capture_kit_vendor'

    ## useless fields
    # 'cases.demographic.ethnicity', # less useful than race; 'not hispanic or latino'
    # 'platform', # generic Illumina
    # 'data_type', # 'Aligned Reads'
    # 'annotations.entity_id', # same as cases.case_id
    # 'annotations.case_id', # same as cases.case_id
    # 'msi_score', # empty
    # 'msi_status', # empty
    # 'proportion_reads_mapped', # empty
    # 'proportion_targets_no_coverage' # empty
    # 'cases.lost_to_followup', # empty
    # 'cases.samples.annotations.entity_id', # empty
    # 'cases.follow_ups.molecular_tests.mismatch_repair_mutation', # empty
    # 'cases.follow_ups.molecular_tests.ploidy', # empty
    # 'analysis.metadata.read_groups.target_capture_kit_version', # empty
    # 'analysis.metadata.read_groups.RIN', # empty
    # 'analysis.metadata.read_groups.instrument_model', # empty
    # 'analysis.metadata.read_groups.read_group_qcs.per_tile_sequence_quality', # PASS/FAIL/WARN
    # 'analysis.metadata.read_groups.read_group_qcs.per_base_sequence_content', # PASS/FAIL/WARN
    # 'analysis.metadata.read_groups.read_group_qcs.per_base_sequence_quality', # PASS/FAIL/WARN
    # 'analysis.metadata.read_groups.read_group_qcs.per_sequence_gc_content', # PASS/FAIL/WARN
    # 'analysis.metadata.read_groups.read_group_qcs.per_sequence_quality_score', # PASS/FAIL/WARN
    # 'analysis.metadata.read_groups.read_group_qcs.workflow_type', # empty
    # 'analysis.metadata.read_groups.read_group_qcs.workflow_version', # empty
    # 'analysis.metadata.read_groups.read_group_qcs.fastq_name', # empty
    # 'analysis.metadata.read_groups.read_group_qcs.basic_statistics', # empty
    # 'downstream_analyses.output_files.proportion_reads_mapped', # empty
    # 'downstream_analyses.output_files.proportion_targets_no_coverage', # empty
    # 'archive.archive_id', # empty
    # 'proportion_coverage_30X', # empty
    # 'average_base_quality', # empty
    # 'mean_coverage', # empty
    # 'cases.sample_ids', # emtpy
    # 'cases.samples.tissue_type', # typically 'not reported'
    # 'cases.samples.distance_normal_to_tumor', # empty
    # 'cases.samples.biospecimen_anatomic_site', # empty
    # 'cases.exposures.pack_years_smoked', # empty
    # 'cases.exposures.alcohol_history', # typically 'not reported'
]
fields = ','.join(fields)

filters = {
    'op': 'and',
    'content': [
        {
            'op': 'or',
            'content': [
                {
                    'op': 'in',
                    'content': {
                        'field': 'cases.project.project_id',
                        'value': 'TCGA-COAD'
                    }
                },
                {
                    'op': 'in',
                    'content': {
                        'field': 'cases.project.project_id',
                        'value': 'TCGA-READ'
                    }
                }
            ]
        },
        {
            'op': 'in',
            'content': {
                'field': 'cases.primary_site',
                'value': ['Colon', 'Rectum', 'Rectosigmoid junction']
            }
        },
        {
            'op': 'in',
            'content': {
                'field': 'experimental_strategy',
                'value': ['WXS']
            }
        },
        {
            'op': 'in',
            'content': {
                'field': 'data_category',
                'value': ['Sequencing Reads']
            }
        }
    ]
}
filters = json.dumps(filters)

params = {
    'filters': filters,
    'fields': fields,
    'format': 'TSV',
    'size': '1500'
    }

response = requests.get(files_endpt, params = params)
object = io.StringIO(response.content.decode('utf-8'))
files_res = pd.read_table(object)

# request total sequences
fields = [
    'file_id',
    'analysis.metadata.read_groups.read_group_qcs.total_sequences'
]
fields = ','.join(fields)

params = {
    'filters': filters,
    'fields': fields,
    'format': 'TSV',
    'size': '1500'
    }

response = requests.get(files_endpt, params = params)
object = io.StringIO(response.content.decode('utf-8'))
seq_bloks = pd.read_table(object)

cols = list(seq_bloks.columns)
pattern = re.compile(r'analysis.metadata.read_groups.\d*.read_group_qcs.0.total_sequences') # identify all possible seq groups and take first of paired-end seq counts
seq_ttl = seq_bloks[list(filter(pattern.match, cols))].sum(axis=1).astype(int)
f_id = seq_bloks['file_id']
frame = {'file_id': f_id, 'total_seq': seq_ttl}
total_seq = pd.DataFrame(frame)

# merge files_res and total_seq
files_res = files_res.merge(total_seq, on='file_id')

# request exome capture kits
index = files_res['file_id'].tolist()
idx_siz = 150
idx_blk = []
idx = []
for i in range(0, len(index)):
    if i == 0:
        idx.append(index[i])
    elif i % idx_siz == 0:
        idx_blk.append(idx)
        idx = [index[i]]
    elif i == (len(index)-1):
        idx.append(index[i])
        idx_blk.append(idx)
    else:
        idx.append(index[i])

fields = [
    'file_id',
    'analysis.metadata.read_groups.target_capture_kit_name',
    'analysis.metadata.read_groups.target_capture_kit_vendor'
]
fields = ','.join(fields)

kits = pd.DataFrame()

for i in idx_blk:
    filters = {
        'op': 'in',
        'content': {
            'field': 'file_id',
            'value': i
        }
    }
    filters = json.dumps(filters)
    params = {
        'filters': filters,
        'fields': fields,
        'format': 'TSV',
        'size': '150'
    }
    response = requests.get(files_endpt, params = params)
    object = io.StringIO(response.content.decode('utf-8'))
    x = pd.read_table(object)
    select = ['file_id', 'analysis.metadata.read_groups.0.target_capture_kit_name', 'analysis.metadata.read_groups.0.target_capture_kit_vendor']
    x = x[select]
    kits = kits.append(x)

kits.columns = ['file_id', 'capture_kit_name', 'capture_kit_vendor']

# merge files_res and capture kits
files_res = files_res.merge(kits, on='file_id')

# format files_res
column_dict = {
    'cases.0.samples.0.is_ffpe': 'is_ffpe',
    'cases.0.exposures.0.height': 'height',
    'cases.0.tissue_source_site.name': 'hospital',
    'cases.0.project.project_id': 'project_id',
    'cases.0.demographic.year_of_birth': 'birth_year',
    'id': 'uuid',
    'cases.0.samples.0.sample_type': 'tissue_type',
    'cases.0.demographic.days_to_death': 'days_to_death',
    'cases.0.demographic.vital_status': 'vital_status',
    'cases.0.samples.0.sample_id': 'sample_id',
    'cases.0.exposures.0.weight': 'weight',
    'cases.0.demographic.race': 'race',
    'cases.0.demographic.age_at_index': 'age_at_index',
    'cases.0.demographic.gender': 'sex',
    'cases.0.case_id': 'case_id',
    'cases.0.primary_site': 'primary_site',
    'cases.0.exposures.0.bmi': 'bmi'
}

column_order = ['file_id', 'uuid', 'file_name', 'sample_id', 'case_id', 'submitter_id', 'project_id', 'primary_site', 'tissue_type', 'birth_year', 'vital_status', 'age_at_index', 'days_to_death', 'height', 'weight', 'bmi', 'race', 'sex', 'hospital', 'is_ffpe', 'capture_kit_name', 'capture_kit_vendor', 'experimental_strategy', 'data_format', 'data_category', 'total_seq', 'file_size', 'md5sum', 'state']

files_res = files_res.rename(column_dict, axis='columns')
files_res = files_res[column_order]





msi_results <- GDCprepare_clinic(query, "msi")




import io
import json
import os
import pandas as pd
import requests
import re
import subprocess
import glob
import xml.etree.ElementTree as ET


data_home = '/scratch/chd5n/aneuploidy/'
anno_path = data_home + 'raw-data/annotations/'


## get uuids of files with the msi annotations from legacy server
## download said files separately
## parse those files to extract msi annotations for each subject

## get uuids
files_leg_endpt = 'https://api.gdc.cancer.gov/legacy/files/'

fields = [
    'file_name',
    'file_id',
    'md5sum',
    'file_size',
    'state'
]
fields = ','.join(fields)

filters = {
    'op':'and',
    'content':[
        {'op':'or',
        'content':[
            {'op':'in',
            'content':{
                'field':'cases.project.project_id',
                'value':'TCGA-COAD'
                }
            },
            {'op':'in',
            'content':{
                'field':'cases.project.project_id',
                'value':'TCGA-READ'
                }
            }
        ]
        },
        {'op':'and',
        'content':[
            {'op':'in',
            'content':{
                'field':'files.data_category',
                'value':'Other'
                }
            },
            {'op':'in',
            'content':{
                'field':'files.data_type',
                'value':'Auxiliary test' ## 'Microsatellite instability' also exists
                }
            },
            {'op':'in',
            'content':{
                'field':'files.access',
                'value':'open'
                }
            }
        ]
        }
    ]
}
filters = json.dumps(filters)

params = {
    'filters': filters,
    'fields': fields,
    'format': 'TSV',
    'size': '1500'
    }

response = requests.get(files_leg_endpt, params = params)
object = io.StringIO(response.content.decode('utf-8'))
files_res = pd.read_table(object)

# ## create manifest object and write to file
#     ## can do this with files_res alone
#     ## id, filename, md5, size, state (state listed as released here, validated from gdc data portal)
# select = ['file_id', 'file_name', 'md5sum', 'file_size', 'state']
# manifest = files_res[select]
# manifest.columns = ['id', 'filename', 'md5', 'size', 'state']
# manifest = manifest.sort_values(by=['id'])
# os.chdir(anno_path)
# manifest.to_csv('msi_manifest.tsv', sep='\t', index=False)


## download
home = os.environ['HOME']
gdc_home = home + '/Apps/gdc-client/'
aneuploidy_home = '/scratch/chd5n/aneuploidy/'
anno_home = aneuploidy_home + 'raw-data/annotations/'
seq_home = aneuploidy_home + 'raw-data/sequencing/'
dest = anno_home + 'msi/'

file_count = 0
for uuid in files_res.id[0:10]:
    cmd = ' '.join([gdc_home + 'gdc-client download',uuid,'-d',dest])
    subprocess.call(cmd, shell=True)
    file_count = file_count + 1
    print(' '.join([uuid,'downloaded']))

print(' '.join([str(file_count),'files downloaded']))

# pat = 'auxiliary:mononucleotide_and_dinucleotide_marker_panel_analysis_status'
#
for uuid in files_res.id[0:10]:
    path = dest + uuid + '/*.xml'
    fn = glob.glob(path)[0]
#
#     with open(fn,'r') as in_f:
#         for in_line in in_f:
#             if re.search(pat, in_line) != None:
#
#                 print(re.search(r'[>a-zA-Z<]', in_line))


for uuid in files_res.id[0:10]:
    path = dest + uuid + '/*.xml'
    fn = glob.glob(path)[0]
    tree = ET.parse(fn)
    root = tree.getroot()
    msi_status = root[1][4][0][1].text
    subject_id = root[1][0].text
    print('\t'.join([subject_id,msi_status]))


bcr patient uuid
root[1][3].text

bcr patient barcode
root[1][0].text


for child in root:
    print(child.tag, child.attrib)

for neighbor in root.iter('neighbor'):
    print(neighbor.attrib)

root = ET.fromstring(country_data_as_string)



nationwidechildrens.org_auxiliary.TCGA-DM-A1D8.xml

https://api.gdc.cancer.gov/data/

xpath <- "//auxiliary:microsatellite_instability_test_result"

grep 'auxiliary:mononucleotide_and_dinucleotide_marker_panel_analysis_status'
