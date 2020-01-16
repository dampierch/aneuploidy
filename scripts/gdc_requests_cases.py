#!/usr/bin/env python

    # this script queries the gdc archive via the search and retrieve api and
    # returns the cases_res object (results from cases endpoint query)


import io
import json
import pandas as pd
import requests


def set_fields():
    '''
    set fields for extraction from endpoint
    '''
    fields = [
        'submitter_id',  ## good
        'case_id',  ## good
        'disease_type',  ## good
        'diagnoses.tumor_stage'  ## good
        # 'sample_ids' ## ok

        # "submitter_analyte_ids", #?
        # "analyte_ids", #?
        # "portion_ids", #?
        # "submitter_portion_ids", #?

        # 'diagnoses.tumor_grade',  ## 'not reported'
        # 'exposures.alcohol_history',  ## 'Not Reported'
        # 'exposures.alcohol_intensity',  ## empty
        # 'exposures.cigarettes_per_day',  ## empty
        # 'exposures.years_smoked'  ## empty
        # "days_to_index", # empty
        # 'samples.biospecimen_anatomic_site', # empty
        # 'samples.distance_normal_to_tumor', # empty
        # 'samples.annotations.case_id', # empty
    ]
    fields = ','.join(fields)
    return fields


def set_filters():
    '''
    set filters to target relevant units on endpoint
    '''
    filters = {
        'op':'and',
        'content':[
            {'op':'or',
            'content':[
                {'op':'in',
                'content':{
                    'field':'project.project_id',
                    'value':'TCGA-COAD'
                }
                },
                {'op':'in',
                'content':{
                    'field':'project.project_id',
                    'value':'TCGA-READ'
                }
                }
            ]
            },
            {'op':'in',
            'content':{
                'field':'primary_site',
                'value':['Colon', 'Rectum', 'Rectosigmoid junction']
            }
            },
            {'op':'in',
            'content':{
                'field':'files.experimental_strategy',
                'value':'WXS'
            }
            },
            {'op':'in',
            'content':{
                'field':'files.data_category',
                'value':'Sequencing Reads'
            }
            }
        ]
    }
    filters = json.dumps(filters)
    return filters


def set_params(filters,fields):
    '''
    set parameters for https get request to endpoint
    '''
    params = {
        'filters': filters,
        'fields': fields,
        'format': 'TSV',
        'size': '1000'
    }
    return params


def get_results(endpoint,params):
    '''
    given an endpoint and parameters, execute https GET request and build
    cases_res dataframe with subject info
    '''
    response = requests.get(endpoint, params=params)
    object = io.StringIO(response.content.decode('utf-8'))
    results = pd.read_table(object)
    return results


def main():
    '''
    get all case_id's from TCGA-COAD and READ
        > 630 unique case_id entries w/o experimental_strategy and data_category filters
        > 608 unique case_id entries (and 1303 files) using filters
    '''
    endpoint = 'https://api.gdc.cancer.gov/cases'
    filters = set_filters()
    fields = set_fields()
    params = set_params(filters,fields)
    cases_res = get_results(endpoint,params)
    cases_res.columns = ['subject_id', 'stage', 'case_id', 'id', 'disease_type']  ## change submitter_id to subject_id for consistency; change diagnoses.0.tumor_stage to stage
    select = ['subject_id', 'case_id', 'disease_type', 'stage']
    cases_res = cases_res[select]
    return cases_res


cases_res = main()
