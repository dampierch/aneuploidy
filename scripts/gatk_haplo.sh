#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50000
#SBATCH --time=6:00:00
#SBATCH --partition=standard
#SBATCH --account=chd5n_alloc
#SBATCH --mail-type=END
#SBATCH --mail-user=chd5n@virginia.edu

## wrp uses #SBATCH --cpus-per-task 16
## HaplotypeCaller appears to look for 4 threads
## runs relatively quickly even with 1 thread

## intersection took about 25 min, <1 GB
## variant call took 24 min to 3h, ~2GB on 1 thread
## process took >10GB on 4 threads for select samples
## process on 15GB slow for samples that went quickly on 10GB (e.g. AF 3400)
## should perhaps choose a memory limit at least as large as largest normal exome bam file (36GB set 1)
## with 50GB, no errors, max 12.3GB, 3h
## outlier: TCGA-AG-4001 normal is 90GB and took 5h (maxrss still 12GB)


module load gcc/7.1.0 bedtools/2.26.0 samtools/1.9 gatk/4.0.0.0

REF_FILE=/scratch/chd5n/Reference_genome/GRCh38.d1.vd1.fa
DATA_DIR=/scratch/chd5n/aneuploidy/raw-data/sequencing/crunch/
exome_bed_file=/scratch/chd5n/aneuploidy/exome-kits/SureSelect_CRE_V2_hg38/S30409818_Regions.bed
raw_bam=$1

cd ${DATA_DIR}

base_name=${raw_bam%.*}
echo `date` :: ${base_name} start
echo haplotype :: ${base_name}
exome_bam=${base_name}_normal_exome.bam

if [ ! -e ${exome_bam} ]; then
    bedtools intersect -a ${raw_bam} -b ${exome_bed_file} > ${exome_bam} ## intersection results in exome bam 60% of original
    samtools index ${exome_bam}
    echo `date` :: ${base_name} intersect done
fi

java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=1 -jar /apps/software/standard/core/gatk/4.0.0.0/gatk-package-4.0.0.0-local.jar HaplotypeCaller -R ${REF_FILE} --input ${exome_bam} -L ${exome_bed_file} -O ${base_name}.snp.indel.vcf
echo `date` :: ${base_name} done

echo "seff ${SLURM_JOBID}"
echo "sacct -o reqmem,maxrss,elapsed,alloccpus -j ${SLURM_JOBID}"
