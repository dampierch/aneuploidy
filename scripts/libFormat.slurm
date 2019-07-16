#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=30000
#SBATCH --time=5:00:00
#SBATCH --partition=standard
#SBATCH --account=cphg_caseylab
#SBATCH --mail-type=ALL
#SBATCH --mail-user=chd5n@virginia.edu


# This script finds the number of paired-end reads in a given BAM file and
# returns a dichotomous variable (single vs paired) to indicate the sequencing
# library format for the sample. If there are >0 paired-end reads, the sample is
# paired, and if there are 0 paired-end reads, the sample is single. This is
# useful so that we indicate the correct library format to functions that act
# on the files.

# Process takes about 6 mins per BAM file for WXS.


##################
# generic set up #
##################

module load samtools

home=/home/chd5n/
project=aneuploidy
scripts=${home}${project}/scripts/
data_master=/scratch/chd5n/${project}/
data_raw=${data_master}raw-data/
seq=${data_raw}sequencing/

SATID=${SLURM_ARRAY_TASK_ID}
subd=$(( 10 * ( ${SATID} - 1 ) ))
subq=$(( 10 * ${SATID} ))

################
# index set up #
################

if [ ${SATID} -eq 1 ]; then
  find ${seq} -mindepth 2 -maxdepth 2 -name *.bam | sed -e "${subq}q" > ${data_raw}format_idx_${SATID}.temp
elif [ ${SATID} -gt 1 ]; then
  find ${seq} -mindepth 2 -maxdepth 2 -name *.bam | sed -e "1,${subd}d;${subq}q" > ${data_raw}format_idx_${SATID}.temp
fi

###########
# execute #
###########

> ${seq}libFormat_${SATID}.tsv
for n in $(<${data_raw}format_idx_${SATID}.temp); do
  d=`echo ${n} | awk -v FS='/' '{print $(NF-1)}'`
  f=`echo ${n} | awk -v FS='/' '{print $NF}' | awk -v FS='_' '{print $1}'`
  c=`samtools view -c -f 1 ${n}`
  if [[ ${c} -eq 0 ]]; then
    lib=single
  elif [[ ${c} -gt 0 ]]; then
    lib=paired
  else
    echo "error with ${n}"
  fi
  printf "${d}\t${f}\t${lib}\n" >> ${seq}libFormat_${SATID}.tsv
done
rm ${data_raw}format_idx_${SATID}.temp
