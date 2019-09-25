#!/bin/bash

  # this script downloads files from gdc with gdc-client using file uuid

# first, we try 10 pairs
gdc_home=~/Apps/gdc-client/
aneuploidy_home=/scratch/chd5n/aneuploidy/
anno_home=${aneuploidy_home}raw-data/annotations/
seq_home=${aneuploidy_home}raw-data/sequencing/
dest=${seq_home}first_10_pairs/
token=${seq_home}gdc-user-token.2019-09-25T22_54_34.618Z.txt
file_info=${anno_home}coad-read.file_info

echo just doing first 20 samples today
printf "start download `date`\n"

file_count=0
for n in `cut -f 4 ${file_info} | head -n 20`; do
    ${gdc_home}gdc-client download ${n} -t ${token} -d ${dest}
    file_count=$(( ${file_count} + 1 ))
    echo ${n} downloaded
done
echo ${file_count} files downloaded

printf "stop download `date`\n"
