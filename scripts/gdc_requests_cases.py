#!/usr/bin/env python

    # this script queries the gdc archive via the search and retrieve api and
    # returns the cases_res object (results from cases endpoint query)

import io
import json
import os
import pandas as pd
import requests

data_home = '/scratch/chd5n/aneuploidy/'
anno_path = data_home + 'raw-data/annotations/'

# get all case_id's from TCGA-COAD and READ
    ## finds 630 unique case_id entries w/o experimental_strategy and data_category filters
    ## finds 608 unique case_id entries (and 1303 files) using filters
cases_endpt = 'https://api.gdc.cancer.gov/cases'

fields = [
    'submitter_id', # good
    'case_id', # good
    'disease_type', # good
    # 'sample_ids' # ok

    # "submitter_analyte_ids", #?
    # "analyte_ids", #?
    # "portion_ids", #?
    # "submitter_portion_ids", #?

    # "days_to_index", # empty
    # 'samples.biospecimen_anatomic_site', # empty
    # 'samples.distance_normal_to_tumor', # empty
    # 'samples.annotations.case_id', # empty
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
                        'field': 'project.project_id',
                        'value': 'TCGA-COAD'
                    }
                },
                {
                    'op': 'in',
                    'content': {
                        'field': 'project.project_id',
                        'value': 'TCGA-READ'
                    }
                }
            ]
        },
        {
            'op': 'in',
            'content': {
                'field': 'primary_site',
                'value': ['Colon', 'Rectum', 'Rectosigmoid junction']
            }
        },
        {
            'op': 'in',
            'content': {
                'field': 'files.experimental_strategy',
                'value': ['WXS']
            }
        },
        {
            'op': 'in',
            'content': {
                'field': 'files.data_category',
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
    'size': '1000'
    }

response = requests.get(cases_endpt, params = params)
object = io.StringIO(response.content.decode('utf-8'))
cases_res = pd.read_table(object)
