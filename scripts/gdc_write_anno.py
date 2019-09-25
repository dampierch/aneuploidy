#!/usr/bin/env python

    # this script creates manifest and pheno objects from gdc query results
    # and writes them to file for future use in download or parsing; the key
    # field offered by cases endpoint is the subject_id (submitter_id)

import io
import json
import os
import pandas as pd
import requests

data_home = '/scratch/chd5n/aneuploidy/'
anno_path = data_home + 'raw-data/annotations/'

from gdc_requests_files import files_res
from gdc_requests_cases import cases_res

# rename cases_res.submitter_id to subject_id, add to files_res, match on case_id
cases_res.columns = ['subject_id', 'case_id', 'id', 'disease_type']
select = ['subject_id', 'case_id', 'disease_type']
cases_res = cases_res[select]
files_res = files_res.merge(cases_res, on='case_id')

# create manifest object and write to file
    ## can do this with files_res alone
    ## id, filename, md5, size, state (state listed as released here, validated from gdc data portal)
select = ['file_id', 'file_name', 'md5sum', 'file_size', 'state']
manifest = files_res[select]
manifest.columns = ['id', 'filename', 'md5', 'size', 'state']
manifest = manifest.sort_values(by=['id'])
os.chdir(anno_path)
manifest.to_csv('manifest.tsv', sep='\t', index=False)

# create pheno annotations object and write to file
select = ['file_id', 'file_name', 'sample_id', 'case_id', 'submitter_id', 'subject_id', 'project_id', 'disease_type', 'primary_site', 'tissue_type', 'birth_year', 'vital_status', 'age_at_index', 'days_to_death', 'height', 'weight', 'bmi', 'race', 'sex', 'hospital', 'is_ffpe', 'total_seq', 'capture_kit_name', 'capture_kit_vendor']
pheno = files_res[select]
pheno = pheno.sort_values(by=['file_id'])
os.chdir(anno_path)
pheno.to_csv('pheno.tsv', sep='\t', index=False)
