# How to use Broad's GATK resources

## Prepare FASTA file as reference

### dict
* Picard tools [CreateSequenceDictionary function guide](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.6.0/picard_sam_CreateSequenceDictionary.php)
* both .fasta and .fasta.gz are supported (gzip)
* also consider this [advisory on Picard syntax transition](https://github.com/broadinstitute/picard/wiki/Command-Line-Syntax-Transition-For-Users-(Pre-Transition))
* process took 1 min, 4.4 GB for TCGA's reference genome

```
# old syntax
java -jar picard.jar CreateSequenceDictionary \
      R=reference.fasta \
      O=reference.dict

# new syntax
java -jar picard.jar CreateSequenceDictionary -R GRCh38.d1.vd1.fa -O GRCh38.d1.vd1.dict

# on Rivanna
module load picard/2.20.6
java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary R=GRCh38.d1.vd1.fa O=GRCh38.d1.vd1.dict
## OR
java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary -R GRCh38.d1.vd1.fa -O GRCh38.d1.vd1.dict
```


### fai
* Samtools [faidx function guide](http://www.htslib.org/doc/samtools-1.2.html)
* input file can be compressed in the BGZF format
* process took < 1 min for TCGA's reference genome

```
# syntax
samtools faidx Homo_sapiens_assembly18.fasta

# on rivanna
module load samtools/1.9
samtools faidx GRCh38.d1.vd1.fa
```

## Call germline SNPs and indels

### HaplotypeCaller
* [documentation page](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php)
