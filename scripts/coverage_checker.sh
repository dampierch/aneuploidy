#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6000
#SBATCH --time=48:00:00
#SBATCH --partition=standard
#SBATCH --account=chd5n_alloc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=chd5n@virginia.edu


  ## this script should be run from:
  ## run_coverage_checker.py with *.file_set specified
  ##
  ## this script takes three command line arguments from run_coverage_checker.py
  ## (subject_id, crunch_path + file_name_norm, crunch_path + file_name_tum)
  ## and finds coverage over a given reference exome for each bam file
  ##
  ## coverage_checker.py appears to use a LOT of time but minimal memory
  ## long-term goal is to parallelize, but for now will use first set results
  ## A6-5661 took ~37h 990MB, AF-3400 took ~15h 660MB


module load anaconda/5.2.0-py3.6

script_dir=~/projects/aneuploidy/scripts/
crunch_dir=/scratch/chd5n/aneuploidy/raw-data/sequencing/crunch/
refex_bed=/scratch/chd5n/aneuploidy/exome-kits/SureSelect_CRE_V2_hg38/S30409818_Regions.bed
subject_id=${1}
normal_pre=${2%.*}
tumor_pre=${3%.*}

echo "################################################################################"
echo
echo starting `date`
echo $1 `echo $2 | awk -v FS="/" '{print $NF}'` `echo $3 | awk -v FS="/" '{print $NF}'`

## get coverage
echo start :: coverage_checker.py :: `date`
${script_dir}coverage_checker.py --bam ${2} \
  --refex_bed ${refex_bed} \
  > ${crunch_dir}${subject_id}_normal_coverage.bed
${script_dir}coverage_checker.py --bam ${3} \
  --refex_bed ${refex_bed} \
  > ${crunch_dir}${subject_id}_tumor_coverage.bed
echo done :: coverage_checker.py :: `date`

echo `date` :: ${base_name} done
echo "seff ${SLURM_JOBID}"
echo "sacct -o reqmem,maxrss,elapsed,alloccpus -j ${SLURM_JOBID}"
echo
echo "################################################################################"
echo
