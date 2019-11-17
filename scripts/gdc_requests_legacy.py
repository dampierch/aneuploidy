#!/usr/bin/env python


    ## this script queries the gdc legacy archive via the search and retrieve
    ## api and returns the msi_status object (from files endpoint on legacy)

    ## 1. get uuids of xml files with the msi annotations from legacy server
    ## 2. download each file
    ## 3. parse files to extract msi annotations for each subject


import io
import json
import os
import pandas as pd
import requests
import re
import subprocess
import glob
import xml.etree.ElementTree as ET


## global variables
home = os.environ['HOME']
gdc_home = home + '/Apps/gdc-client/'
aneuploidy_home = '/scratch/chd5n/aneuploidy/'
anno_home = aneuploidy_home + 'raw-data/annotations/'
seq_home = aneuploidy_home + 'raw-data/sequencing/'


def set_filters():
    '''
    set filters for gdc legacy files endpoint search; these filters determine
    which type of data is requested and retrieved; json format
    '''
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


def get_results(endpoint,filters):
    '''
    request xml file_id entities for msi results from TCGA-COAD and READ
    '''
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
    response = requests.get(endpoint, params=params)
    object = io.StringIO(response.content.decode('utf-8'))
    results = pd.read_table(object)
    return results


def download_xml_uuid(files_res,gdc_home,dest):
    '''
    download xml files one at a time by uuid
    '''
    file_count = 0
    for uuid in files_res.id:
        cmd = ' '.join([gdc_home + 'gdc-client download',uuid,'-d',dest])
        subprocess.call(cmd, shell=True)  ## download xml
        print(' '.join([uuid,'downloaded']))
        file_count = file_count + 1
    print(' '.join([str(file_count),'files downloaded']))


def download_xml_manifest(files_res,gdc_home,dest):
    '''
    1. create manifest object; 2. write to file; 3. use manifest for download
    '''
    select = ['file_id', 'file_name', 'md5sum', 'file_size', 'state']
    manifest = files_res[select]
    manifest.columns = ['id', 'filename', 'md5', 'size', 'state']
    manifest = manifest.sort_values(by=['id'])
    out_file = anno_home + 'msi_manifest.tsv'
    manifest.to_csv(out_file, sep='\t', index=False)
    cmd = ' '.join([gdc_home + 'gdc-client download','-m',out_file,'-d',dest])
    subprocess.call(cmd, shell=True)
    print('manifest downloaded')


def parse_xml(files_res,dest):
    '''
    parse xml files to extract msi status
    '''
    fields = ['subject_id','msi_status']
    msi_dict = {}
    msi_dict['subject_id'] = []
    msi_dict['msi_res'] = []
    file_count = 0
    for uuid in files_res.id:
        pattern = dest + uuid + '/*.xml'
        fn = glob.glob(pattern)[0]
        tree = ET.parse(fn)  ## parse xml
        root = tree.getroot()
        subject_id = root[1][0].text
        msi_status = root[1][4][0][1].text
        msi_dict['subject_id'].append(subject_id)
        msi_dict['msi_res'].append(msi_status)
        file_count = file_count + 1
    print(' '.join([str(file_count),'files parsed']))
    msi_res = pd.DataFrame.from_dict(msi_dict)
    return msi_res


def main():
    endpoint = 'https://api.gdc.cancer.gov/legacy/files/'
    filters = set_filters()
    files_res = get_results(endpoint,filters)
    dest = anno_home + 'msi/'
    download_xml_manifest(files_res,gdc_home,dest)
    msi_res = parse_xml(files_res,dest)
    return msi_res

## bcr patient uuid
# root[1][3].text
## bcr patient barcode
# root[1][0].text

msi_res = main()
