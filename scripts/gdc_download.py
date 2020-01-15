#!/usr/bin/env python3


  ## this script downloads files from gdc with gdc-client using file uuid;
  ## executes download on front end


import subprocess
from datetime import datetime
import os
import argparse
import glob


def get_uuids(in_file):
    '''
    read all uuids into list
    '''
    with open(in_file,'r') as in_f:
        uuids = []
        for in_line in in_f:
            uuids.append(in_line.strip('\n').split('\t')[-1])
    return uuids


def parse_uuids(uuids,set_size):
    '''
    break uuids into manageable sets for download
    '''
    uuid_sets = []
    uuid_set = []
    for i in range(0, len(uuids)):
        if i == 0:
            uuid_set.append(uuids[i])
        elif i % set_size == 0:
            uuid_sets.append(uuid_set)
            uuid_set = [uuids[i]]
        elif i == (len(uuids)-1):
            uuid_set.append(uuids[i])
            uuid_sets.append(uuid_set)
        else:
            uuid_set.append(uuids[i])
    return uuid_sets


def download_files(uuid_sets,set_num,gdc_home,dest,token,uuids_undone):
    '''
    download files specified in uuid_set on frontend
    '''
    print(' '.join(['download attempt start',str(datetime.now())]))
    file_count = 0
    if set_num == 12:
        targets = uuids_undone
    elif set_num == 13:
        targets = ['97401027-5a1a-498f-b256-004b28b9db44','afd08843-c655-45d4-873d-a31c826b1727']
    else:
        targets = uuid_sets[set_num]
    for uuid in targets:
        cmd = ' '.join([gdc_home + 'gdc-client download',uuid,'-d',dest,'-t',token])
        subprocess.call(cmd, shell=True)
        file_count = file_count + 1
        print(uuid + ' attempted')
    print(str(file_count) + ' files attempted')
    print(' '.join(['download attempt end',str(datetime.now())]))


def write_downloaded(downloaded_file,uuid_sets,set_num,uuids_undone):
    '''
    append to record the uuid set downloaded
    '''
    with open(downloaded_file,'w+') as out_file:
        line1 = ' '.join(['most recent set number downloaded:',str(set_num)]) + '\n'
        if set_num == 12:
            line2 = ' '.join(uuids_undone) + '\n'
        elif set_num == 13:
            line2 = ' '.join(['97401027-5a1a-498f-b256-004b28b9db44','afd08843-c655-45d4-873d-a31c826b1727']) + '\n'
        else:
            line2 = ' '.join(uuid_sets[set_num]) + '\n'
        out_file.write(line1+line2)


def main(set_num):
    '''
    parse uuids and download selected set
    '''
    ## set variables
    home = os.environ['HOME']
    gdc_home = home + '/Apps/gdc-client/'
    aneuploidy_home = '/scratch/chd5n/aneuploidy/'
    anno_home = aneuploidy_home + 'raw-data/annotations/'
    seq_home = aneuploidy_home + 'raw-data/sequencing/'
    in_file = anno_home + 'coad-read.file_info'
    dest = seq_home + 'current_set/'
    token = seq_home + 'gdc-user-token.2020-01-06T02_18_41.663Z.txt'
    set_size = 100
    set_num = int(set_num)
    downloaded_file = seq_home + 'latest_download.txt'
    ## take actions
    uuids = get_uuids(in_file)
    uuid_sets = parse_uuids(uuids,set_size)
    paths_done = glob.glob(dest+'*/')
    uuids_done = [i.split('/')[-2] for i in paths_done]
    uuids_undone = [uuid for uuid in uuids if uuid not in uuids_done]
    download_files(uuid_sets,set_num,gdc_home,dest,token,uuids_undone)
    write_downloaded(downloaded_file,uuid_sets,set_num,uuids_undone)


parser = argparse.ArgumentParser(description='Download set of UUIDs')
parser.add_argument('--set_num',help='uuid set number',action='store',dest='set_num')
args = parser.parse_args()
main(args.set_num)


## testing to see which sets are complete
# set_size = 100
# uuid_sets = parse_uuids(uuids,set_size)
#
# set_dict = {}
# converse = {}
# for i in range(0,12,1):
#     if i in set_dict:
#         pass
#     else:
#         set_dict[i] = []
#         converse[i] = []
#
# for i in range(0,12,1):
#     set_dict[i] = [uuid for uuid in uuids_done if uuid in uuid_sets[i]]
#     converse[i] = [uuid for uuid in uuids_done if uuid not in uuid_sets[i]]
#
# len_ids = [len(set_dict[key]) for key in set_dict]
# [80, 20, 100, 100, 100, 100, 100, 100, 100, 100, 100, 74]
# len_ids = [len(converse[key]) for key in converse]
# [994, 1054, 974, 974, 974, 974, 974, 974, 974, 974, 974, 1000]
