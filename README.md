# Aneuploidy in CRC

An investigation into the drivers of chromosomal instability (CIN) in colorectal cancer (CRC) using The Cancer Genome Atlas (TCGA) whole exome sequencing (WXS) data from tumor samples and matched non-tumor samples.

## Pipeline
* For the live pipeline, see [Snakefile](scripts/Snakefile)
* For comments on development, see [Pipeline Readme](readme_pipeline.md)

## Analysis v0.1
* DGE using all samples, most vs least aneuploid
* DGE using MSS samples, most vs least aneuploid

## Analysis v0.2
* add normal SD curve to aneuploid rank plot to help with interpretation - DONE
* VAAST is proprietary - ABANDON - NEED TO FIND ALTERNATIVE
* gene expression correlations with MSS and MSI-L - DONE
* share het-sites directory with Pankaj - DONE
  1. `/scratch/chd5n/aneuploidy/`
  2. `/scratch/chd5n/aneuploidy/hetsites-data/r-cnts/`

## Analysis v0.3
* classify chromosome arms using 3 state HMM
* find patterns

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

| set | download | assemble | variants | hetsites | plots | store | LOH | trouble |
| :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: | :--: |
| 1 | x | x | x | x | x | x |  |  |
| 2 | x | x | x | x | x | x |  |  |
| 3 | x | x | x | x | x | x |  |  |
| 4 | x | x | x | x | x | x |  |  |
| 5 | x | x | x | x | x | x |  |  |
| 2v2 | x | x | x | x | x | x |  |  |
| 3v2 | x | x | x | x | x | x |  |  |
| 4v2 | x | x | x | x | x | x |  | TCGA-AG-4005 |
| 5v2 | x | x | x | x | x | x | TCGA-CK-5912, TCGA-DM-A1DA | TCGA-AG-4001 |
| 6v2 | x | x | x | x | x | x |  |  |
| 7v2 | x | x | x | x | x | x | TCGA-DY-A1DG | TCGA-AG-3999 |
| 8v2 | x | x | x | x | x | x |  | TCGA-AG-4007 |
| 9v2 | x | x | x | x | x | x | TCGA-DM-A1HA | TCGA-AG-A00H | ... problem with TCGA-AA-3866
| 10v2 | x | x | x | x | x | x |  |  |
| 11v2 | x | x | x | x | x | x |  |  |
| 12v2 | x | x | x | x | x | x |  |  |

Total: 587 pairs (1174 WXS files) processed

## Background
* Pfister et al., 2018: [Identification of Drivers of Aneuploidy in Breast Tumors](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5997284/)
1. Ried et al. (review), 2019: [The landscape of genomic copy number alterations in colorectal cancer and their consequences on gene expression levels and disease outcome](https://www.sciencedirect.com/science/article/pii/S0098299719300354)
2. Taylor et al., 2018: [Genomic and Functional Approaches to Understanding Cancer Aneuploidy](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6028190/)
3. Sansregret & Swanton (review), 2017: [The Role of Aneuploidy in Cancer Evolution](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5204330/)
4. Andor et al., 2016: [Pan-cancer analysis of the extent and consequences of intra-tumor heterogeneity](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4830693/)

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
