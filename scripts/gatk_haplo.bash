#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=128000
#SBATCH --time=8:00:00
#SBATCH --partition=standard
#SBATCH --account=chd5n_alloc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=chd5n@virginia.edu

## wrp also uses #SBATCH --cpus-per-task 16


module load gcc/7.1.0 bedtools/2.26.0 samtools/1.9 gatk/4.0.0.0

REF_FILE=/scratch/chd5n/Reference_genome/GRCh38.d1.vd1.fa
DATA_DIR=/scratch/chd5n/aneuploidy/raw-data/sequencing/crunch/
exome_bed_file=/scratch/chd5n/aneuploidy/exome-kits/SureSelect_CRE_V2_hg38/S30409818_Regions.bed
raw_bam=$1

cd ${DATA_DIR}

# #idx=`cut -f 2 $1`
# idx=`cut -f 2 /scratch/chd5n/aneuploidy/raw-data/annotations/coad-read_2019-09-26.file_set`
# for raw_bam in ${idx}; do

base_name=${raw_bam%.*}
echo `date` :: ${base_name} start
echo haplotype :: ${base_name}
exome_bam=${base_name}_normal_exome.bam

if [ ! -e ${exome_bam} ]; then
    bedtools intersect -a ${raw_bam} -b ${exome_bed_file} > ${exome_bam}
    samtools index ${exome_bam}
    echo `date` :: ${base_name} intersect done
fi

java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=1 -jar /apps/software/standard/core/gatk/4.0.0.0/gatk-package-4.0.0.0-local.jar HaplotypeCaller -R ${REF_FILE} --input ${exome_bam} -L ${exome_bed_file} -O ${base_name}.snp.indel.vcf_L
echo `date` :: ${base_name} done
