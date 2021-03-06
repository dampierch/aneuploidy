#!/usr/bin/env python3

## unzip_de.py :: unzip gzip-compressed count matrix files downloaded from gdc
    ## usage: python unzip_de.py --set_name {rna_set_name}


import subprocess
from datetime import datetime
import os
import argparse


def get_uuids(in_file):
    '''
    read uuids and file names into lists
    '''
    with open(in_file,'r') as in_f:
        uuids = []
        fnames = []
        header = True
        for in_line in in_f:
            if header:
                header = False
                continue
            uuids.append(in_line.strip('\n').split('\t')[-3])
            fnames.append(in_line.strip('\n').split('\t')[-4])
    return uuids,fnames


def unzip_files(uuids,fnames,count_home):
    '''
    check uuid and file name concordance and unzip files
    '''
    print(' '.join(['unzip start',str(datetime.now())]))
    new_file_count = 0
    old_file_count = 0
    for uuid in uuids:
        if os.path.exists('%s/%s' % (count_home + uuid, fnames[new_file_count][:-3])):
            old_file_count = old_file_count + 1
            new_file_count = new_file_count + 1
            continue
        elif os.path.exists('%s/%s' % (count_home + uuid, fnames[new_file_count])):
            cmd = ' '.join(['gunzip','%s/%s' % (count_home + uuid, fnames[new_file_count])])
            subprocess.call(cmd, shell=True)
            # print(cmd)  ## for testing
            print(uuid + ' attempted')
        else:
            print('%s %s' % ('file name discrepancy for',uuid))
        new_file_count = new_file_count + 1
    print(str(old_file_count) + ' files already decompressed')
    print(str(new_file_count) + ' files considered')
    print(' '.join(['unzip end',str(datetime.now())]))


def main(set_name):
    '''
    parse uuids and unzip files
    '''
    ## set variables
    home = os.environ['HOME']
    aneuploidy_home = '/scratch/chd5n/aneuploidy/'
    anno_home = aneuploidy_home + 'raw-data/annotations/'
    count_home = aneuploidy_home + 'raw-data/counts/'
    in_file = anno_home + 'rna_set_' + set_name + '.tsv'
    ## take actions
    uuids,fnames = get_uuids(in_file)
    unzip_files(uuids,fnames,count_home)


parser = argparse.ArgumentParser(description='Gunzip set of UUIDs')
parser.add_argument('--set_name',help='uuid set name',action='store',dest='set_name')
args = parser.parse_args()
main(args.set_name)
