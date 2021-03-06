#!/usr/bin/env python

    ## this script creates manifest and pheno objects from gdc query results
    ## and writes them to file for future use in download or parsing; the key
    ## field offered by cases endpoint is the subject_id (submitter_id)


import io
import json
import os
import pandas as pd


from gdc_requests_files import files_res
from gdc_requests_cases import cases_res
from gdc_requests_legacy import msi_res


def merge_results(files_res,cases_res,msi_res):
    '''
    add cases_res to files_res matching on case_id (this adds subject_id);
    then add msi_res to cases_res matching on subject_id
    '''
    results = files_res.merge(cases_res, on='case_id')   ## inner join by default
    results = results.merge(msi_res, how='left', on='subject_id')   ## left join
    return results


def make_manifest(gdc_anno,anno_home):
    '''
    create manifest object from gdc_anno and write to file;
    manifest fields: id, filename, md5, size, state
    '''
    select = ['file_id', 'file_name', 'md5sum', 'file_size', 'state']
    manifest = gdc_anno[select]
    manifest.columns = ['id', 'filename', 'md5', 'size', 'state']
    manifest = manifest.sort_values(by=['id'])
    manifest.to_csv(anno_home + 'manifest.tsv', sep='\t', index=False)


def make_pheno(gdc_anno,anno_home):
    '''
    create pheno annotations object from gdc_anno and write to file
    '''
    select = ['file_id', 'file_name', 'sample_id', 'case_id', 'submitter_id', 'subject_id', 'project_id', 'disease_type', 'primary_site', 'tissue_type', 'stage', 'msi_status', 'birth_year', 'vital_status', 'age_at_index', 'days_to_death', 'height', 'weight', 'bmi', 'race', 'sex', 'hospital', 'is_ffpe', 'total_seq', 'capture_kit_name', 'capture_kit_vendor', 'file_name_rna', 'file_id_rna']
    pheno = gdc_anno[select]
    pheno = pheno.sort_values(by=['file_id'])
    pheno.to_csv(anno_home + 'pheno.tsv', sep='\t', index=False)


def main():
    '''
    make a manifest and phenotype annotation for dataset of interest
    '''
    anno_home = '/scratch/chd5n/aneuploidy/raw-data/annotations/'
    gdc_anno = merge_results(files_res,cases_res,msi_res)
    make_manifest(gdc_anno,anno_home)
    make_pheno(gdc_anno,anno_home)


main()
