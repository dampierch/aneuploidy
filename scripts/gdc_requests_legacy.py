#!/usr/bin/env python


    ## this script queries the gdc legacy archive via the search and retrieve
    ## api and returns the msi_status object (from files endpoint on legacy)


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
endpoint = 'https://api.gdc.cancer.gov/legacy/files/'



def set_filters():
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
    return filters

fields = [
    'file_name',
    'file_id',
    'md5sum',
    'file_size',
    'state'
]
fields = ','.join(fields)
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
