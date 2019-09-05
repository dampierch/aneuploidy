import requests
import json
import os

# get all case_id's from TCGA-COAD and READ
    ## this gives 630 unique case_id entries
    ## expecting 608 based on website, maybe some have no data
    ## we do get 608 unique case_id entries (and 1303 files) using files endpoint
cases_endpt = 'https://api.gdc.cancer.gov/cases'

fields = [
    'submitter_id',
    'case_id',
    'primary_site',
    'disease_type',
    'project.project_id'
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

with open('/scratch/chd5n/aneuploidy/test_request.txt', 'w+') as f:
    f.write(response.content.decode('utf-8'))
