# Aneuploidy in CRC

An investigation into the drivers of CIN in CRC using TCGA WXS.

## Pipeline

### WRP::request file info from GDC
* uses JSON request to files endpoint via curl
* identifies files by filter for cases.submitter_id (subject_id)

#### CHD::request file info from GDC
* uses JSON request to files, cases, and legacy endpoints via python
* identifies files by filter for cases.project.project_id, cases.primary_site, experimental_strategy, and data_category
* generates a manifest for bulk download with gdc-client and a pheno file for parsing
* legacy endpoint gives msi_status; cases endpoint gives submitter_id and disease_type; files endpoint gives rest
* [gdc_requests_cases.py](scripts/gdc_requests_cases.py)
* [gdc_requests_files.py](scripts/gdc_requests_files.py)
* [gdc_requests_legacy.py](scripts/gdc_requests_legacy.py)
* [gdc_write_anno.py](scripts/gdc_write_anno.py)

### WRP::parse_tcga_info.py hall_tcga_t10b10.tab > hall_tcga_t10b10.file_info
* generates file with following fields: sample_id | sample_type | file_name | file_id
* for normal/tumor pairs (in that order)
* unclear how it deals with duplicated samples

#### CHD::modify parse_tcga_info.py to take as input pheno.tsv and generate samples.file_info
* i do so in the following script and of duplicates i take max sequences
* output is coad-read.file_info
* [gdc_parse_info.py](scripts/gdc_parse_info.py)

### WRP::nohup ~/ncbi/gdc_download_file_info.sh ~/ncbi/hall_tcga_t61-80.file_info > gdc_down_t61-80.log 2> gdc_down_t61-80.err &
* bash script feeds file uuid from parse*.py output to gdc-client instead of manifest
* this allows more precise download but does not appear to permit the ever-ellusive sbatch download

#### CHD::use WRP download strategy (can't use bulk download anyway due to not enough storage space)
* first try 9/25/2019 on first ten pairs (20 samples) in coad-read.file_info
* old: `gdc_download.bash`
* update: [gdc_download.py](scripts/gdc_download.py)

```
dt=`date +"%Y-%m-%d"`
nohup bash ~/projects/aneuploidy/scripts/gdc_download.bash > ~/projects/aneuploidy/logs/gdc_download_${dt}.out 2> ~/projects/aneuploidy/logs/gdc_download_${dt}.err &
```

* seems to have worked

### WRP::move_gdc_tcga_files.py hall_tcga_t10b10 hall_tcga_t10b10.file_info > hall_tcga_t10b10.file_set 2> hall_tcga_t10b10.file_errors
* `*.file_set` contains: subject_id, normal_file.bam, tumor_file.bam
* need to move bai files as well, because later pysam will need them

#### CHD::check and assemble files for analysis
* old: `assemble_gdc_files.py`
* [gdc_assemble_files.py](scripts/gdc_assemble_files.py)

### WRP:Identify heterozygous sites and count reads there
* de-duplication is assumed to have been done by GDC
* requires a `*.file_set` file to identify samples
* takes normal_file.bam and runs gatk_haplo.sh via sbatch

```
for n in `cut -f 2 hall_tcga_t9b10.file_set`; do sbatch run_gatk_haplo.sh $n ; done
```

* run_gatk_haplo.sh does three things:
  1. bedtools intersect with exome.bed
  2. samtools index
  3. java gatk.jar HaplotypeCaller
* the script takes three input files:
  1. reference.fa file (UCSC hg38; requires dict and index)
  2. reference exome (from manufacturer, with liftOver from hg19 to hg38)
  3. bam file (from .file_set)
* the script generates two output files:
  1. `*_normal_ed.bam`
  2. `*.snp.indel.vcf_L`

```
start_paired_hets.py hall_tcga_t61.file_set
```

* calls find_count_hets_tumor_pair_gdc.sh, which in turn calls find_hetsites.py, count_het_freqs2.py, het_cnts2R.py
* find_hetsites.py does three things:
  1. reads normal sample vcf
  2. selects variant positions with AF=0.5 and good depth
  3. writes bed file with position of baseline heterozygous sites
  4. output: `*_het.bed`
  5. fields: `'CHROM', 'POS0', 'POS', 'name', 'QUAL', 'strand', 'info_str'`
* count_het_freqs2.py does three things:
  1. reads an unfiltered bam file (tumor or normal)
  2. reads bed file `*_het.bed` with position of baseline (i.e. normal) heterozygous sites
  3. writes vcf-like output with allele counts in bam file (tumor or normal) at each location in bed file
  4. output: `*.het_cnts2`
  5. the allele counts for the normal bam may be redundant w/r/t the `*.snp.indel.vcf_L` allele counts...we will see
* het_cnts2R.py does two things:
  1. reads allele counts at heterozygous locations in tumor and normal `*.het_cnts2`
  2. prepares allele count information for easy input into R
  3. output: `out_fields = ['sample', 'chrom', 'pos', 't_type', 'ref', 'ref_cnt', 'alt', 'alt_cnt', 'fract']`
* het_cnts2R_n.py does a better job of checking to ensure that the same chromosome positions are being merged

#### CHD::gatk_haplo
* which exome.bed to use: either broad, generic or sample-specific; going to try broad, generic first
* i downloaded SureSelect Clinical Research Exome V2 (67.3 Mb, hg38) from agilent catalogue...no liftOver required!
* will use S30409818_Regions.bed as exome reference

```
/scratch/chd5n/aneuploidy/exome-kits/SureSelect_CRE_V2_hg38/S30409818_Regions.bed
```

* which reference genome.fa to use: either generic latest hg38 from ucsc or GDC-specific; going to try GDC-specific first

```
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/chromFa.tar.gz .
rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz .
```

* [Heng Li recommends](http://lh3.github.io/2017/11/13/which-human-reference-genome-to-use):
* `ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz`
* [GDC used](https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files) a combo of Heng Li's rec plus decoys and viruses

```
# for Li's
rsync -avzhP rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz /scratch/chd5n/Reference_genome/ ## took 9 sec

# for GDC (use gdc-client)
gdc_home=~/Apps/gdc-client/
dest=/scratch/chd5n/Reference_genome/
${gdc_home}gdc-client download 254f697d-310d-4d7d-a27b-27fbf767a834 -d ${dest}  ## took 7 sec
cd ${dest}254f697d-310d-4d7d-a27b-27fbf767a834
tar xvzf GRCh38.d1.vd1.fa.tar.gz
# now resides at following path:
/scratch/chd5n/Reference_genome/GRCh38.d1.vd1.fa
```

* create dict and index for reference genome
* run gatk_haplo
  1. [run_gatk_haplo.py](scripts/run_gatk_haplo.py)
  2. [gatk_haplo.sh](scripts/gatk_haplo.sh)

```
cd ~/projects/aneuploidy/scripts
file_set=/scratch/chd5n/aneuploidy/raw-data/annotations/coad-read_2019-09-26.file_set
for file in `cut -f 2 ${file_set}`; do sbatch gatk_haplo.bash ${file}; done
```

* output is vcf_L (TCGA-A6-5661-10A-01D-1650-10_Illumina_gdc_realn is large outlier)

```
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  TCGA-T9-A92H-10A-01D-A370-10
chr1    17385   .       G       A       176.77  .       AC=1;AF=0.500;AN=2;BaseQRankSum=-2.255;ClippingRankSum=0.000;DP=39;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=39.04;MQRankSum=-1.732;QD=4.53;ReadPosRankSum=0.201;SOR=0.099       GT:AD:DP:GQ:PL  0/1:30,9:39:99:205,0,1033
chr1    69511   .       A       G       3492.77 .       AC=2;AF=1.00;AN=2;DP=164;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=35.95;QD=21.83;SOR=1.565       GT:AD:DP:GQ:PL  1/1:0,160:160:99:3521,463,0

## INFO README
AC=allele count in genotypes, for each alt allele
AF=allele freq, for each alt allele (theoretical)
AN=total number of alleles (always 2 in this set)
DP=approx read depth, some reads filtered
DS=any downsampling?
QD=QualByDepth
MLEAC=max likelihood expectation for allele count, for each alt allele
MLEAF=max likelihood expectation for allele freq, for each alt allele

## FORMAT README
GT=genotype 0/1
AD=allelic depths for ref/alt alleles 30,9
DP=approx read depth 39 (sum of allelic depths)
GQ=genotype quality 99
PL=normalized, phred-scaled likelihoods for genotypes 205,0,1033
```

* [1000G reference page for genotype format in vcf](https://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40/)
* AF is theoretical based on called genotype, so always 0.5 or 1; empirical genotype can be calculated from AD field
* homozygous genotypes are sometimes called even with reads supporting heterozygosity; presumably based on read and alignment quality

```
cut -f 8 TCGA-T9-A92H-10A-01D-A370-10_Illumina_gdc_realn.snp.indel.vcf_L | awk -v FS=";" 'NR>2807 {print $2}' | sort | uniq -c | less
22058 AF=0.500
   72 AF=0.500,0.500
10580 AF=1.00

cut -f 8 TCGA-A6-2680-11A-01D-1554-10_Illumina_gdc_realn.snp.indel.vcf_L | awk -v FS=";" 'NR>2807 {print $2}' | sort | uniq -c | less
22664 AF=0.500
   44 AF=0.500,0.500
11304 AF=1.00

cut -f 8 TCGA-A6-2680-11A-01D-1554-10_Illumina_gdc_realn.snp.indel.vcf_L | awk -v FS=";" 'NR>2807 {print $3}' | sort | uniq -c | less
awk -v FS="\t" 'NR>2807 {print $0}' TCGA-A6-2680-11A-01D-1554-10_Illumina_gdc_realn.snp.indel.vcf_L | less
grep 'AF=1.00' TCGA-A6-2680-11A-01D-1554-10_Illumina_gdc_realn.snp.indel.vcf_L | cut -f 10 | awk -v FS=":" '{print $2}' | sort | uniq -c | less
```

#### CHD::allele counts
* start_paired_hets.py `*.file_set` :: see [run_het_counter.py](scripts/run_het_counter.py)
* find_count_hets_tumor_pair_gdc.sh :: see [het_counter.sh](scripts/het_counter.sh)
* find_hetsites.py :: see [find_hetsites.py](scripts/find_hetsites.py); includes modification to accommodate genotypes without reference allele as well as QD threshold instead of depth and QUAL
* count_het_freqs2.py :: see [count_hetalleles.py](scripts/count_hetalleles.py); includes modification for more intuitive loop over pileupcolumns as well as cov_thresh
* het_cnts2R.py :: see [hetcnts_2R.py](scripts/hetcnts_2R.py); includes minor modification with get_hetcnts_list function
* must remember to make scripts executable with `#!/usr/bin/env python3` at top and `chmod +x` at unix command line
* with dp_thresh in count_hetalleles set to x, first_10_pairs had:

```
wc -l *_normal_errcnts.bed
wc -l *_normal_hetcnts.bed
grep 'normal missing' *_R_missing.err | wc -l
wc -l *_tumor_errcnts.bed
wc -l *_tumor_hetcnts.bed
grep 'tumor missing' *_R_missing.err | wc -l
```

| dp_thresh | tissue type | errcnts | hetcnts | R missing |
| :--: | :--: | :--: | :--: | :--: |
| 100 | normal | 13760 | 69617 | 4282 |
| 100 | tumor | 36830 | 46547 | 28574 |
| 0 | normal | 37 | 83340 | 10 |
| 0 | tumor | 2130 | 81247 | 2103 |
| 50 | normal | 43 | 83334 | 17 |
| 50 | tumor | 14101 | 69276 | 14041 |

* with qd_thresh in find_hetsites set to 5 and cov_thresh in count_hetalleles set to 20, first_10_pairs had:

| tissue type | errcnts | hetcnts | R missing |
| :--: | :--: | :--: | :--: |
| normal | 15764 | 155493 | 3503 |
| tumor | 28820 | 142437 | 29298 |

##### CHD:optimal dp_thresh (i.e. cov_thresh)
* will be related to coverage; see [run_coverage_checker.py](scripts/run_coverage_checker.py), [coverage_checker.sh](scripts/coverage_checker.sh), and [coverage_checker.py](scripts/coverage_checker.py)
* these scripts work but take too long...need to parallelize the depth counter...
* code from this [BioStars post](https://www.biostars.org/p/275974/#276179) should be helpful...
* more importantly, GATK has DiagnoseTargets; see [coverage_analysis.py](scripts/coverage_analysis.py)


### wrp::plots
* these scripts require an 'hg38.chr1-XY.sizes_cumm' file with cumulative chromosome lengths
* `Rscript --vanilla plot_haplo_gg_chrom2.R $n` -- plot two panel histogram (density) and one-panel chromosome plot  -- try this first
* `plot_karyo_gg_list.R`   -- plot using a file with a list of TCGA id's, makes karyoploteR plot
* `plot_karyo_gg.R`        -- plot an individual TCGA id, makes karyoploteR plot


#### CHD::plots
* plot_haplo_gg_chrom2r.R :: see [run_plotter.py](scripts/run_plotter.py) and [plot_haplo_chrom.R](scripts/plot_haplo_chrom.R)
* karyoploteR...


### smallest file in test set
/scratch/chd5n/aneuploidy/raw-data/sequencing/crunch/TCGA-AF-3400-01A-01D-1989-10_gapfillers_Illumina_gdc_realn.bam

### msi status
* solved with [gdc_requests_legacy.py](scripts/gdc_requests_legacy.py)

### snakemake pipeline
* should be ready to go: [Snakefile](scripts/Snakefile)...
* trying set 1
  1. start 2019-11-20: download worked, 20 in 2 hours >> fe
  2. assemble_files worked (3-4sec) >> fe
  3. call_variants worked (1h 40m) >> sbatch
  4. count_hetsites worked (50m) >> sbatch
  5. make_density_plots worked after bug fix (9-12sec) >> fe
  6. store_hetsite_data working after bug fix (20sec) >> fe

### sets processed

| set | download | assemble | variants | hetsites | plots | store |
| :--: | :--: | :--: | :--: | :--: | :--: | :--: |
| 1 | x | x | x | x | x | x |
| 2 | x | x | x | x | x | x |
| 3 | x | x | x | x | x | x |
| 4 | x | x | x | x | x | x |
| 5 | x | x | x | x | x | x |
| 2v2 | x | x | x | x | x | x |
| 3v2 | x | x | x | . |  |  |
| 4v2 |  |  |  |  |  |  |
| 5v2 |  |  |  |  |  |  |
| 6v2 |  |  |  |  |  |  |
| 7v2 |  |  |  |  |  |  |
| 8v2 |  |  |  |  |  |  |
| 9v2 |  |  |  |  |  |  |
| 10v2 |  |  |  |  |  |  |
| 11v2 |  |  |  |  |  |  |
| 12v2 |  |  |  |  |  |  |


## early thoughts

### Sample selection
1. From the [GDC data portal](https://portal.gdc.cancer.gov/), we select all colon, rectosigmoid, and rectum samples from TCGA-COAD and TCGA-READ without a filter for tumor type (i.e. we do not care whether the diagnosis was adenocarcinoma or cystic neoplasm).
2. We add all WXS files to our cart, and we find there are a total of 1,303 files from 608 unique cases (i.e. subjects); this data takes up 31.31 TB of disk space (including index files).
3. From our cart, we download the following annotation files:
  * biospecimen.cart.2019-07-13.json
  * biospecimen.cart.2019-07-13.tar.gz
  * clinical.cart.2019-07-13.json
  * clinical.cart.2019-07-13.tar.gz
  * gdc_sample_sheet.2019-07-13.tsv
  * metadata.cart.2019-07-13.json
4. Using gdc_sample_sheet.2019-07-13.tsv along with [samples.py](scripts/samples.py), we explore the basic characteristics of our cohort.
5. We can mostly replicate this dataset with [gdc_requests_files.py](scripts/gdc_requests_files.py).
  * We can identify the sequencing files but can't easily capture the annotation files, although most of the important fields can be captured with [gdc_write_anno.py](scripts/gdc_write_anno.py)

| Sample type | *n* |
| :--: | :--: |
| Primary Tumor | 628 |
| Blood Derived Normal | 560 |
| Solid Tissue Normal | 112 |
| Recurrent Tumor | 2 |
| Metastatic | 1 |

6. Among identifiers in the sample sheet, there are 1,259 unique 'Sample ID's, 608 unique 'Case ID's, 1,303 unique 'File Name's and 'File ID's.

#### Normal tissue
1. We start by examining the normal tissue from which we will call SNVs.
2. There are 672 samples of normal tissue WXS.
3. There are 601 unique cases among those samples; 71 samples are non-unique.
4. There are 68 cases with multiple tissue samples, 65 with just two and three with a total of three samples each, accounting for all 71 non-unique samples (139 files).
5. We may have a soft preference for blood-derived normal WXS over adjacent normal tissue, but the choice at this point is arbitrary between duplicate blood-derived normals due to lack of quality parameters.
6. There are 533 cases with unique tissue samples, 477 blood-derived and 56 solid tissue (i.e. adjacent normal).

#### Tumor tissue
1. We next examine the tumor tissue we will test for functional aneuploidy.
2. There are 628 samples of tumor tissue WXS.
3. There are 594 unique cases among those samples; 34 samples are non-unique.
4. There are 28 cases with multiple tissue samples, 22 with just two and six with a total of three samples each, accounting for all 34 non-unique samples (62 files).
5. At this point, the choice is arbitrary among the duplicates, as all are Primary Tumor, and we do not have quality parameters.
6. There are 566 cases with unique tissue samples, all primary tumor.

#### Pilot
1. We create a pilot set of five cases from the set of unique normal samples and match them with their five tumor counterparts. We use [samples.py](scripts/samples.py) to generate a pilot manifest.
2. We use the [gdc-client](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool) to download BAM files to a `/scratch` directory.
3. After checking library format with [libFormat.slurm](scripts/libFormat.slurm), we use [Biobambam2](https://www.sanger.ac.uk/science/tools/biobambam) to revert the BAM files to FASTQ and [pigz](https://zlib.net/pigz/) to compress them.

### Variants

1. From our normal WXS data, we want to call germline SNVs. We could also call somatic mutations, but the TCGA VCF files should suffice for that.
2. We start with our subset of five normal tissue samples that are unique and have a unique tumor sample partner. There are 510 of these unique pairs in total.
3. Our [download script](download.slurm) does not seem to work on the compute nodes, so we run it on a login node. We expect the download to be about 240 GB (average file size of [31 210 000 000 000 bytes / 1303 files] = 24 GB * 10 files). The actual file size is 218 GB for 10 BAM files.
4. We call variants with GATK-HC

### Variant calling pipelines
* [Kumaran, 2019](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2928-9) finds BWA or Novoalign with DeepVariant or SAMTools give best SNV results (GIAB gold standard)
* [Hwang, 2015](https://www.nature.com/articles/srep17875?report=reader) finds BWA-MEM with SAMTools gives best SNV results; also Samtools tends to add reference alleles and thereby overcall heterozygous SNVs (GIAB gold standard)
* [Pirooznia, 2014](https://humgenomics.biomedcentral.com/articles/10.1186/1479-7364-8-14), using BWA with realignment/recalibration, finds GATK-UG outperforms SAMTools mpileup and GATK-HC better than GATK-UG (Sanger gold standard)... may be due to GATK outperformance for indels
* [Liu, 2013](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0075619), using BWA, finds GATK outperforms SAMTools (Sanger gold standard)... may be due to GATK outperformance for indels

### Genome in a Bottle (NIST)
* [Zook, 2014](https://www.nature.com/articles/nbt.2835) introduces GIAB (uses GATK-HC + GATK-UG + Cortex for gold standard variant calls)

### Exome capture kits
* [Wang, 2018](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0204912) reports that exome capture kits bias results for limited number of genes; [Github](https://github.com/TheJacksonLaboratory/GDCSlicing) may have code to extract exome capture kit from metadata
* [WouterDeCoster](https://www.biostars.org/p/220939/) suggests using primary target coordinates for analysis (as opposed to capture targets, i.e. the bait); this is consistent with target_capture_kit_target_region url provided in TCGA metadata from GDC
* [Baylor kit](https://sequencing.roche.com/en/products-solutions/by-category/target-enrichment/hybridization/seqcap-ez-hgsc-vcrome.html)

| Kit | *n* |
| :--: | :--: |
| SeqCap EZ HGSC VCRome | 878 |
| SureSelect Human All Exon 38 Mb v2 | 186 |
| Gapfiller_7m | 94 |
| SeqCap EZ Human Exome Library v2.0 | 49 |
| VCRome V2.1 | 47 |
| VCRomeV2.1-PKv1 | 31 |
| Custom V2 Exome Bait, 48 RXN X 16 tubes | 2 |
| SeqCap EZ Exome V2.0 | 1 |
| NaN | 15 |


# Scratch space

## R code for GDC TCGA API
        # if(json$data$pagination$count == 0) {
        #     url <- getGDCquery(project = proj,
        #                        data.category = data.category,
        #                        data.type = data.type,
        #                        legacy = legacy,
        #                        workflow.type = NA,
        #                        platform = NA,
        #                        file.type = file.type,
        #                        experimental.strategy = experimental.strategy,
        #                        files.access = access,
        #                        sample.type = sample.type)
        #     json  <- tryCatch(
        #         getURL(url,fromJSON,timeout(600),simplifyDataFrame = TRUE),
        #         error = function(e) {
        #             message(paste("Error: ", e, sep = " "))
        #             message("We will retry to access GDC!")
        #             fromJSON(content(getURL(url,GET,timeout(600)), as = "text", encoding = "UTF-8"), simplifyDataFrame = TRUE)
        #         }
        #     )
        # }
