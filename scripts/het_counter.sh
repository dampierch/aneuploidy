#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2000
#SBATCH --time=2:00:00
#SBATCH --partition=standard
#SBATCH --account=chd5n_alloc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=chd5n@virginia.edu

################
## this script should be run from:
## run_het_counter.py with *.file_set specified
##
## this script takes three command line arguments from run_het_counter.py
## (subject_id, crunch_path + file_name_norm, crunch_path + file_name_tum)
## and does three things: (1) finds het sites from normal vcf (2) captures AF
## at each het site found in (1) for normal and tumor (3) prepares AF data for
## downstream processing in R
##
## find_hetsites and count_hetalleles appear to use minimal time/mem
## first 10 samples took 22m - 1h; max ~55 MB
## hetcnts_2R.py uses no detectable memory, takes 1 sec per sample


module load anaconda/5.2.0-py3.6

script_dir=~/projects/aneuploidy/scripts/
crunch_dir=/scratch/chd5n/aneuploidy/raw-data/sequencing/crunch/
vcf_suf="snp.indel.vcf"
subject_id=${1}
normal_pre=${2%.*}
tumor_pre=${3%.*}

echo starting `date`
echo $1 `echo $2 | awk -v FS="/" '{print $NF}'` `echo $3 | awk -v FS="/" '{print $NF}'`

## find het sites
echo start finding hetsites `date`
${script_dir}find_hetsites.py ${normal_pre}.${vcf_suf} > ${normal_pre}_hetsites.bed
echo done finding hetsites `date`

## get emperical allele counts
echo start counting alleles `date`
${script_dir}count_hetalleles.py --bam ${2} \
  --hetsites_bed ${normal_pre}_hetsites.bed \
  > ${crunch_dir}${subject_id}_normal_hetcnts.bed \
  2> ${crunch_dir}${subject_id}_normal_errcnts.bed
${script_dir}count_hetalleles.py --bam ${3} \
  --hetsites_bed ${normal_pre}_hetsites.bed \
  > ${crunch_dir}${subject_id}_tumor_hetcnts.bed \
  2> ${crunch_dir}${subject_id}_tumor_errcnts.bed
echo done counting alleles `date`

## prepare data for R
echo start preparing R data `date`
${script_dir}hetcnts_2R.py --subject ${subject_id} \
  --normal_hetcnts_bed ${crunch_dir}${subject_id}_normal_hetcnts.bed \
  --tumor_hetcnts_bed ${crunch_dir}${subject_id}_tumor_hetcnts.bed \
  2> ${crunch_dir}${subject_id}_R_missing.err
echo done preparing R data `date`

echo `date` :: ${base_name} done
echo "seff ${SLURM_JOBID}"
echo "sacct -o reqmem,maxrss,elapsed,alloccpus -j ${SLURM_JOBID}"
