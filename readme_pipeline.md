# Bioinformatic Pipeline README

## Query files available for study

To determine what files are available for analysis, a list of files meeting certain criteria was requested from the Genomic Data Commons.

### Original WRP solution
1. JSON request to `files` endpoint via `curl`
2. Identify files by filter for `cases.submitter_id` (subject_id)
* *This method works when you know which samples you want to download*

### CHD solution
1. JSON request to `files`, `cases`, and `legacy` endpoints via Python scripts
2. Identify files by filter for `cases.project.project_id`, `cases.primary_site`, `experimental_strategy`, and `data_category`
3. Generate a manifest for bulk download with `gdc-client` and a pheno file for parsing
* *`legacy` endpoint gives msi_status; `cases` endpoint gives submitter_id and disease_type; `files` endpoint gives all other information*
4. The Python scripts (write_anno is the master script):
* [gdc_requests_cases.py](scripts/gdc_requests_cases.py)
* [gdc_requests_files.py](scripts/gdc_requests_files.py)
* [gdc_requests_legacy.py](scripts/gdc_requests_legacy.py)
* [gdc_write_anno.py](scripts/gdc_write_anno.py)

## Parse and organize available files

To organize the files available, samples were grouped by subject into tumor-normal pairs.

### Original WRP solution
1. Generate file with following fields: sample_id | sample_type | file_name | file_id
2. List normal/tumor pairs (in that order, i.e. normal first)
* *Unclear how to manage duplicated samples (i.e. multiple tumors or normals for single subject)*
3. Executed with the following code:

```
parse_tcga_info.py hall_tcga_t10b10.tab > hall_tcga_t10b10.file_info
```

### CHD solution
1. Modify parse_tcga_info.py to take as input pheno.tsv and generate samples.file_info
2. For duplicates, take sample with max sequences (i.e. reads)
3. Output is `coad-read.file_info`
4. The Python script:
* [gdc_parse_info.py](scripts/gdc_parse_info.py)

## Download files

To download the exome sequencing reads, files of interest were requested by uuid and transferred from Genomic Data Commons servers to local servers using `gdc-client` software.

### Original WRP solution
1. Bash script feeds file uuid from parse_tcga_info.py output to `gdc-client` (does not use manifest)
* *UUID approach allows more precise download (i.e. per file) but does not appear to permit the ever-ellusive `sbatch` download via compute node*

```
nohup ~/ncbi/gdc_download_file_info.sh ~/ncbi/hall_tcga_t61-80.file_info > gdc_down_t61-80.log 2> gdc_down_t61-80.err &
```

### CHD solution
1. Download by UUID instead of manifest (not enough storage space for bulk download)
2. v0.1 2019-09-25 using [SLURM](scripts/download.slurm) with manifest failed
3. v0.2 2019-09-25 on first ten pairs (20 samples) in coad-read.file_info using `gdc_download.bash` worked

```
dt=`date +"%Y-%m-%d"`
nohup bash ~/projects/aneuploidy/scripts/gdc_download.bash > ~/projects/aneuploidy/logs/gdc_download_${dt}.out 2> ~/projects/aneuploidy/logs/gdc_download_${dt}.err &
```

4. v1.0 2019-11-20 Python script created to download batches of UUIDs: [gdc_download.py](scripts/gdc_download.py)

## Check download and assemble files

To verify exome sequencing reads were successfully transferred, BAM files from each pair of tumor-normal samples were checked for existence. For complete pairs, BAM files and associated indices (i.e. BAI files) were assembled in a single directory to facilitate analysis.

### Original WRP solution
1. `*.file_set` contains: subject_id, normal_file.bam, tumor_file.bam
2. Must keep BAI files with BAM files because `pysam` will need them

```
move_gdc_tcga_files.py hall_tcga_t10b10 hall_tcga_t10b10.file_info > hall_tcga_t10b10.file_set 2> hall_tcga_t10b10.file_errors
```

### CHD solution
1. v0.1 check and then assemble files using the now deprecated `assemble_gdc_files.py`
2. v1.0 check and assemble files with [gdc_assemble_files.py](scripts/gdc_assemble_files.py)

## Identify SNPs

To identify heterozygous sites to use for estimations of ploidy, we sought sites of common variation (i.e. single nucleotide polymorphisms) by calling germline variants in bulk exome sequencing reads from normal samples for each subject separately. To select for SNPs likely to have the best coverage, and to save computation time, we called variants in the subset of exome sequencing reads that aligned to a reference exome.

### Original WRP solution
1. De-duplication of DNA sequencing reads is assumed to have been completed by TCGA prior to upload of BAM files to GDC
2. A `*.file_set` file is required to identify samples for processing
3. A normal BAM is selected and variants are called with the GATK HaplotypeCaller algorithm using SLURM
4. Prior to variant calling, input BAM files are filtered to call variants on just the subset of reads that align to a reference exome (given in BED format)
* *Whole exome sequencing is supposed to select DNA templates of exons alone, but off-target sequences are captured as well. Non-exonic regions are biologically suitable for variant calling, but we exclude them because we expect lower coverage of these regions and concomitant lower accuracy in haplotype estimation. By excluding non-exonic regions, we save processing time and exclude SNPs likely to be called with lower confidence.*
5. Filtered BAM file must be indexed prior to calling variants
6. Command:

```
for n in `cut -f 2 hall_tcga_t9b10.file_set`; do sbatch run_gatk_haplo.sh $n ; done
```

#### `run_gatk_haplo.sh` does three things

1. bedtools intersect exome.bed (bam with reference exome)
2. samtools index (filtered bam)
3. java gatk.jar HaplotypeCaller (filtered bam)

#### `run_gatk_haplo.sh` requires three inputs

1. reference genome against which to call variants; FASTA format (UCSC hg38; requires dict and index)
2. reference exome to filter reads and limit haplotype search space; BED format (from manufacturer, with liftOver from hg19 to hg38)
3. aligned reads on which to call variants; BAM format (from .file_set)

#### `run_gatk_haplo.sh` generates two outputs

1. `*_normal_filtered.bam`
2. `*.snp.indel.vcf`

### CHD solution
1. Use same approach: filter reads for intersection with exons and then call variants
2. Choose an appropriate reference exome: SureSelect Clinical Research Exome V2
3. Choose an appropriate reference genome (and create dict and index per [GATK guidelines](readme_gatk.md)): GRCh38.d1.vd1.fa
4. Run bedtools intersect, samtools index, gatk HaplotypeCaller with the following scripts:
* [run_gatk_haplo.py](scripts/run_gatk_haplo.py)
* [gatk_haplo.sh](scripts/gatk_haplo.sh)

#### Choice of reference exome
1. Could use broad, generic exome or a sample-specific exome to match specific exome capture kit used for each sample (see [Exome Capture Kit README](readme_exome_capture.md))
2. The same broad, generic reference exome would be easier to implement and might decrease sample-specific or capture kit-specific biases
3. A reasonable choice is the **SureSelect Clinical Research Exome V2** (67.3 Mb, hg38) from Agilent
4. SureSelect CRE V2 has more total coverage than two most common kits used in COAD and READ, although I am not certain the TCGA kits are subsets of CRE V2
* *SeqCap EZ HGSC VCRome (45.1 Mb, Roche) and SureSelect Human All Exon v2 (38 Mb, Agilent) were the two most common TCGA kits, both probably designed based on hg19*
5. SureSelect CRE V2 has advantage of being in hg38 coordinates, so no liftOver required to use with GDC BAM files and reference genome, which are both GRCh38
6. Chose S30409818_Regions.bed as exome reference (obtained from Agilent)

```
/scratch/chd5n/aneuploidy/exome-kits/SureSelect_CRE_V2_hg38/S30409818_Regions.bed
```

#### Choice of reference genome
1. Could use generic GRCh38 from ENSEMBL or modified semi-generic builds from UCSC, Broad, ENCODE, NCBI
2. GDC has specific version of GRCh38 used for latest TCGA alignments; will use this one
3. Heng Li [recommends](http://lh3.github.io/2017/11/13/which-human-reference-genome-to-use) the so-called "no alt analysis set"

```
ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
```

4. [GDC used](https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files) a combo of Heng Li's recommended version plus decoys and viruses
5. Code for downloading reference genomes:

```
cd /scratch/chd5n/downloads/

## ucsc
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/chromFa.tar.gz .
rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz .

## no alt analysis set (9 sec download)
rsync -avzhP rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz .

## gdc via gdc-client (7 sec download)
gdc_home=~/Apps/gdc-client/
dest=/scratch/chd5n/downloads/
${gdc_home}gdc-client download 254f697d-310d-4d7d-a27b-27fbf767a834 -d ${dest}
cd ${dest}254f697d-310d-4d7d-a27b-27fbf767a834
tar xvzf GRCh38.d1.vd1.fa.tar.gz
```

#### Notes from first try
1. Command used:

```
cd ~/projects/aneuploidy/scripts
file_set=/scratch/chd5n/aneuploidy/raw-data/annotations/coad-read_2019-09-26.file_set
for file in `cut -f 2 ${file_set}`; do sbatch gatk_haplo.bash ${file}; done
```

2. Sample output: (note TCGA-A6-5661-10A-01D-1650-10_Illumina_gdc_realn is an outlier)

```
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  TCGA-T9-A92H-10A-01D-A370-10
chr1    17385   .       G       A       176.77  .       AC=1;AF=0.500;AN=2;BaseQRankSum=-2.255;ClippingRankSum=0.000;DP=39;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=39.04;MQRankSum=-1.732;QD=4.53;ReadPosRankSum=0.201;SOR=0.099       GT:AD:DP:GQ:PL  0/1:30,9:39:99:205,0,1033
chr1    69511   .       A       G       3492.77 .       AC=2;AF=1.00;AN=2;DP=164;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=35.95;QD=21.83;SOR=1.565       GT:AD:DP:GQ:PL  1/1:0,160:160:99:3521,463,0
```

3. More sample output:

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

4. VCF codes

```
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

* [1000G reference page for genotype format in VCF](https://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40/)
* AF is theoretical based on called genotype, so always 0.5 or 1; empirical genotype can be calculated from AD field
* homozygous genotypes are sometimes called even with reads supporting heterozygosity; presumably based on read and alignment quality

## Identify heterozygous sites

Having identified germline SNPs in normal tissue, we next identified heterozygous sites by selecting the subset of SNPs with sufficient sequencing coverage and quality of variant call for which the estimated allele fraction was 0.5 in normal tissue.

### Original WRP solution
1. Python master script calls shell script that subsequently calls Python scripts to identify heterozygous sites, count them, and prepare them for analysis
2. May be useful to break these steps into two: first the identification and then the counting and prep
3. Master script is `start_paired_hets.py`, shell script is `find_count_hets_tumor_pair_gdc.sh`
4. Heterozygous site identification is performed by a VCF parser script called `find_hetsites.py`
5. Command:

```
start_paired_hets.py hall_tcga_t61.file_set
```

#### `find_hetsites.py` does three things
1. Reads normal sample vcf
2. Selects variant positions with AF=0.5 and good depth (*good depth is unknown, initially set to 100*)
3. Writes bed file with position of baseline heterozygous sites: `*_het.bed`
4. Output fields: `'CHROM', 'POS0', 'POS', 'name', 'QUAL', 'strand', 'info_str'`

### CHD solution
1. Mimic WRP solution with Python master script to call shell script that subsequently calls Python scripts to identify heterozygous sites, count them, and prepare them for analysis
2. May be useful to break these steps into two: first the identification and then the counting and prep
3. Master script is [run_het_counter.py](scripts/run_het_counter.py), shell script is [het_counter.sh](scripts/het_counter.sh)
4. Heterozygous site identification is performed by the modified VCF parser script [find_hetsites.py](scripts/find_hetsites.py)
* *Modified to accommodate genotypes without reference allele (admittedly very few) as well as QD threshold instead of depth and QUAL thresholds*
* *See [here](https://github.com/dampierch/aneuploidy/blob/de710eb8106b0601de5d8f166d2aba9039c28861/scripts/find_hetsites.py#L76) for QD threshold and [here](https://github.com/dampierch/aneuploidy/blob/de710eb8106b0601de5d8f166d2aba9039c28861/scripts/find_hetsites.py#L81) and [here](https://github.com/dampierch/aneuploidy/blob/de710eb8106b0601de5d8f166d2aba9039c28861/scripts/find_hetsites.py#L35) for alt only alleles*

#### Choice of QD threshold
1. QD=QualByDepth value
2. GATK recommends 2 (i.e. exclude variants with QD<2)
3. I use a slightly more conservative threshold of 5 (i.e. exclude variants with QD<5)
4. Optimal threshold for this purpose is unknown; some notes on the subject and some preliminary tests are in [coverage_analysis.py](scripts/coverage_analysis.py)

## Quantify variant allele fractions

Having identified germline heterozygous sites for each subject, we next quantified allele fractions at those sites in normal and tumor tissue. To quantify allele fractions, numbers of unfiltered reads of each allele at all heterozygous sites with adequate coverage were counted. Sites of interest are necessarily bi-allelic in this context because they are polymorphic with respect to the reference genome in normal, diploid cells from individual subjects. Sites of interest thus typically harbor one reference and one alternate allele (for sites that are bi-allelic in the population) but may rarely harbor two different alternate alleles (for sites that are poly-allelic in the population).

### Original WRP solution
1. As noted above, Python master script calls shell script that subsequently calls Python scripts to identify heterozygous sites, count them, and prepare them for analysis
2. Master script is `start_paired_hets.py`, shell script is `find_count_hets_tumor_pair_gdc.sh`
3. Unfiltered reads at heterozygous sites are counted with `pysam` in `count_het_freqs2.py`
4. Read counts are prepared for analysis with `het_cnts2R.py`
5. Command:

```
start_paired_hets.py hall_tcga_t61.file_set
```

#### `count_het_freqs2.py` does three things
1. Reads an unfiltered bam file (tumor or normal)
2. Reads bed file `*_het.bed` with position of germline (i.e. normal) heterozygous sites
3. Counts reads overlapping each allele at each site of interest (from bed) and reports allele counts for sites with sufficient coverage and without too high a percentage of unexpected alleles; writes vcf-like output to `*.het_cnts2`
* *Optimal coverage and expected allele thresholds are unknown.*
* *Allele counts for normal should be equal to AD value in `.snp.indel.vcf_L` but this is untested.*

#### `het_cnts2R.py` does two things
1. Reads allele counts at heterozygous sites from `*.het_cnts2`
2. Prepares allele count information for easy input into R by filtering for ?something? ?sites with paired information (i.e. both normal and tumor counts exist at the site in question)? and by calculating major allele fractions
* *Note `het_cnts2R_n.py` does a better job of checking to ensure that the same chromosome positions are being merged*

### CHD solution
1. As noted above, mimic WRP solution with Python master script to call shell script that subsequently calls Python scripts to identify heterozygous sites, count them, and prepare them for analysis
2. Master script is [run_het_counter.py](scripts/run_het_counter.py), shell script is [het_counter.sh](scripts/het_counter.sh)
3. Unfiltered reads at heterozygous sites are counted with `pysam` in [count_hetalleles.py](scripts/count_hetalleles.py)
* *Modified to use more intuitive loop over pileupcolumns and to include a semi-arbitrary, reasonable, mildly conservative coverage threshold*
4. Read counts are prepared for analysis (i.e. converted to allele fractions) with [hetcnts_2R.py](scripts/hetcnts_2R.py)
* *Modified slightly with `get_hetcnts_list` function*
* *Implementation Note: remember to make scripts executable with `#!/usr/bin/env python3` at top and `chmod +x` at Unix command line*

#### Choice of coverage threshold
1. `dp_thresh` (a.k.a. `cov_thresh`) originally set to 100 in `count_het_freqs2.py` and `het_cnts2R.py`
2. With change to using `qd_thresh` of 5 in [find_hetsites.py](scripts/find_hetsites.py), we can try lowering `cov_thresh` in [count_hetalleles.py](scripts/count_hetalleles.py) to 20
3. Optimal threshold for this purpose is unknown but should depend on empiric coverage; some notes on the subject and some preliminary tests are in [coverage_analysis.py](scripts/coverage_analysis.py) as well as [run_coverage_checker.py](scripts/run_coverage_checker.py), [coverage_checker.sh](scripts/coverage_checker.sh), and [coverage_checker.py](scripts/coverage_checker.py)
* *The actual coverage counting is performed by [coverage_checker.py](scripts/coverage_checker.py)*
* *Implementation Note: these scripts work but take too long; the depth counter should be parallelized*
* *Parallelization of `pysam` is subject of this [Biostars post](https://www.biostars.org/p/275974/#276179)*
* *GATK DepthOfCoverage function or old DiagnoseTargets function is probably what I should use; more in [coverage_analysis.py](scripts/coverage_analysis.py)*

#### Coverage threshold tests
1. Simple commands for preliminary test:

```
wc -l *_normal_errcnts.bed
wc -l *_normal_hetcnts.bed
grep 'normal missing' *_R_missing.err | wc -l
wc -l *_tumor_errcnts.bed
wc -l *_tumor_hetcnts.bed
grep 'tumor missing' *_R_missing.err | wc -l
```

2. With `qual_thresh = 100` and `dp_thresh = 100` in [find_hetsites.py](scripts/find_hetsites.py) and variable `dp_thresh` in [count_hetalleles.py](scripts/count_hetalleles.py), we find the following counts in first_10_pairs set:

| dp_thresh | tissue type | errcnts | hetcnts | R missing |
| :--: | :--: | :--: | :--: | :--: |
| 100 | normal | 13760 | 69617 | 4282 |
| 100 | tumor | 36830 | 46547 | 28574 |
| 0 | normal | 37 | 83340 | 10 |
| 0 | tumor | 2130 | 81247 | 2103 |
| 50 | normal | 43 | 83334 | 17 |
| 50 | tumor | 14101 | 69276 | 14041 |

* Counted heterozygous sites not very sensitive to `dp_thresh` in [count_hetalleles.py](scripts/count_hetalleles.py)

3. With `qd_thresh = 5` in [find_hetsites.py](scripts/find_hetsites.py) and `cov_thresh = 20` in [count_hetalleles.py](scripts/count_hetalleles.py), we find the following counts in first_10_pairs set:

| tissue type | errcnts | hetcnts | R missing |
| :--: | :--: | :--: | :--: |
| normal | 15764 | 155493 | 3503 |
| tumor | 28820 | 142437 | 29298 |

* Counted heterozygous sites were more restricted by original `qual_thresh` and `dp_thresh` in [find_hetsites.py](scripts/find_hetsites.py)
* `qd_thresh = 5` in [find_hetsites.py](scripts/find_hetsites.py) and `cov_thresh = 20` in [count_hetalleles.py](scripts/count_hetalleles.py) is relatively permissive

#### Heterozygous site counts
1. Show distribution of counted heterozygous sites across samples and check correlation with known technical factors (e.g. exome capture kit, source)
2. Also check missing chromosomes problem
3. Use [sum_sites.py](scripts/sum_sites.py) and [sum_sites.R](scripts/sum_sites.R)
4. We see that all samples with fewer than 5000 sites tested are from Gapfiller_7m exome capture kit
5. The maximum sites for a tumor sample with the Gapfiller_7m kit is 2665; the analgous value for normal samples is 2563
6. We see a cluster of samples in the tumor TSV (i.e. Rdata) filtered set of sites with site counts between 5000 and 8000 and about 10 chromosomes included
7. We cannot fix the Gapfiller_7m outliers; we should exclude them
8. We can fix the missing chromosome problem, and we expect this cluster to disappear after fixing the problem with the sorting, intercalating, and filtering script (i.e. hetcnts_2R)

## HMM classifier
* *Under construction* [fit_hmm.py](scripts/fit_hmm.py)

## TO DO
* Fix [hetcnts_2R.py](scripts/hetcnts_2R.py) to make 'phasing' more sensible
* Fix [plot_haplo_chrom.R](scripts/plot_haplo_chrom.R) to simplify syntax
* Finish [fit_hmm.py](scripts/fit_hmm.py)
