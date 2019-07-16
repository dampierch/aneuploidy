#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100000
#SBATCH --time=02:00:00
#SBATCH --partition=standard
#SBATCH --account=cphg_caseylab
#SBATCH	--mail-type=ALL
#SBATCH	--mail-user=chd5n@virginia.edu

# This script converts BAM to FASTQ using Biobambam2. Whereas BEDTools takes
# about 105 min per BAM file and Picard Tools takes about 30 min, Biobambam2
# takes about 5 min for RNASeq.

# Process takes about 10 min per BAM file for WXS.

##################
# generic set up #
##################

home=/home/chd5n/
project=aneuploidy
scripts=${home}${project}/scripts/
data_master=/scratch/chd5n/${project}/
data_raw=${data_master}raw-data/
seq=${data_raw}sequencing/
bbb2_path=${home}Apps/biobambam2/2.0.87-release-20180301132713/x86_64-etch-linux-gnu/bin/
lib_format=${seq}libFormat_pilot.tsv

SATID=${SLURM_ARRAY_TASK_ID}
subd=$(( 10 * ( ${SATID} - 1 ) ))
subq=$(( 10 * ${SATID} ))

################
# index set up #
################

if [ ${SATID} -eq 1 ]; then
  find ${seq} -mindepth 2 -maxdepth 2 -name *.bam | sed -e "${subq}q" > ${data_raw}bbb2_idx_${SATID}.temp
elif [ ${SATID} -gt 1 ]; then
  find ${seq} -mindepth 2 -maxdepth 2 -name *.bam | sed -e "1,${subd}d;${subq}q" > ${data_raw}bbb2_idx_${SATID}.temp
fi

###########
# execute #
###########

for n in $(<${data_raw}bbb2_idx_${SATID}.temp); do
  d=`echo ${n} | awk -v FS='/' '{print $(NF-1)}'`
  f=`echo ${n} | awk -v FS='/' '{print $NF}' | awk -v FS='_' '{print $1}'`
  lib=`awk -v d=${d} -v FS='\t' '$1==d {print $NF}' ${lib_format}`
  if [[ ${lib} == 'single' ]]; then
    ${bbb2_path}/bamtofastq filename=${n} S=${seq}${d}/${f}.fastq
  elif [[ ${lib} == 'paired' ]]; then
    ${bbb2_path}/bamtofastq filename=${n} F=${seq}${d}/${f}_1.fastq F2=${seq}${d}/${f}_2.fastq
  else
    echo "error in ${n}"
  fi
done
rm ${data_raw}bbb2_idx_${SATID}.temp
