#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=30000
#SBATCH --time=10:00:00
#SBATCH --partition=standard
#SBATCH --account=cphg_caseylab
#SBATCH --mail-type=ALL
#SBATCH --mail-user=chd5n@virginia.edu

# This script compresses FASTQ files with pigz, which utilizes multithreading
# to speed up the compression.

# Process took 01:40 hh:mm with pigz for 428 gtex RNASeq fastq files (n=31 in 14 batches)
# Process took 07:36 hh:mm with pigz for 292 sra pub RNASEQ fastq files in series

##################
# generic set up #
##################

home=/home/chd5n/
project=aneuploidy
scripts=${home}${project}/scripts/
data_master=/scratch/chd5n/${project}/
data_raw=${data_master}raw-data/
seq=${data_raw}sequencing/

SATID=${SLURM_ARRAY_TASK_ID}
subd=$(( 20 * ( ${SATID} - 1 ) ))
subq=$(( 20 * ${SATID} ))

################
# index set up #
################

if [ ${SATID} -eq 1 ]; then
  find ${seq} -mindepth 2 -maxdepth 2 -name "*.fastq" | sed -e "${subq}q" > ${data_raw}pigz_idx_${SATID}.temp
elif [ ${SATID} -gt 1 ]; then
  find ${seq} -mindepth 2 -maxdepth 2 -name "*.fastq" | sed -e "1,${subd}d;${subq}q" > ${data_raw}pigz_idx_${SATID}.temp
fi

###########
# execute #
###########

for n in $(<${data_raw}pigz_idx_${SATID}.temp); do
  pigz -9 -p 8 ${n}
done
rm ${data_raw}pigz_idx_${SATID}.temp
exit


# Alternative index for more complex cases:

################
# index set up #
################

batchSize=31
batchNum=14

d=$(( ${batchSize} * ( ${SATID} - 1 ) ))
q=$(( ${batchSize} * ${SATID} ))

#parse full index once
if [ ${SATID} -eq 1 ]; then

  > wait.temp
  printf "task ID = ${SATID}\nwill parse index now...\n"
  find . -type f -name "*.fastq" | sort > idx.temp
  rm wait.temp
  > idxDone.temp
  printf "task ID = ${SATID}\nindex complete, will proceed...\n"

elif [ ${SATID} -gt 1 ]; then

  printf "task ID = ${SATID}\nchecking for index...\n"
  sleep 1
  while [ -e wait.temp ]; do
    printf "task ID = ${SATID}\nwaiting for index...\n"
    sleep 10
  done
  if [ -e idxDone.temp ]; then
    printf "task ID = ${SATID}\nindex complete, will proceed...\n"
  else
    printf "task ID = ${SATID}\nparsing trouble, exiting with exit code 3\n"
    exit 3

  fi
fi

#set batch indices
if [ ${SATID} -eq 1 ]; then
  awk -v FS='/' '{print $NF}' idx.temp | sed -e "${q}q" > gzIdx${SATID}.temp
elif [ ${SATID} -gt 1 ]; then
  awk -v FS='/' '{print $NF}' idx.temp | sed -e "1,${d}d;${q}q" > gzIdx${SATID}.temp
fi

###########
# execute #
###########

for n in $(<gzIdx${SATID}.temp); do
  pigz -9 -p 8 ${n}
done
rm gzIdx${SATID}.temp
