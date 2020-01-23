#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20                       #### set to 20 for BiocParallel
#SBATCH --mem=050000                             #### set to 150 for large BiocParallel
#SBATCH --time=00:30:00                          #### hh:mm:ss
#SBATCH --partition=standard
#SBATCH --account=chd5n_alloc                    #### alternative: cphg_caseylab
#SBATCH --mail-type=END
#SBATCH --mail-user=chd5n@virginia.edu

## variables
script_name=de.R
rna_set_name=$1
sva_pass_num=$2
sva_refine=$3
sva_nsv=$4

## header
pwd; hostname; date
## environments
module purge
# module load gcc/7.1.0 openmpi/3.1.4 R/3.6.0
module load gcc/7.1.0 R/3.5.1
# module load gcc/7.1.0 R/3.6.1
## execution
printf "start ${script_name}\n"
Rscript ${script_name} --args ${rna_set_name} ${sva_pass_num} ${sva_refine} ${sva_nsv}
## footer
printf "end ${script_name}\n"; date
echo "seff ${SLURM_JOBID}"
echo "sacct -o reqmem,maxrss,elapsed,alloccpus,nodelist -j ${SLURM_JOBID}"
