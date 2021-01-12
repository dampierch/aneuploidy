# Bioinformatic Pipeline README

## Query files available for study

To determine what files are available for analysis, a list of files meeting certain criteria was requested from the Genomic Data Commons.

### Original WRP solution
1. JSON request to `files` endpoint via `curl`
2. identify files by filter for `cases.submitter_id` (subject_id)
  * This method works when you know which samples you want to download

### CHD solution
1. JSON request to `files`, `cases`, and `legacy` endpoints via Python scripts
2. identify files by filter for `cases.project.project_id`, `cases.primary_site`, `experimental_strategy`, and `data_category`
  * legacy endpoint gives msi_status; cases endpoint gives submitter_id and disease_type; files endpoint gives rest
  * [gdc_requests_cases.py](scripts/gdc_requests_cases.py)
  * [gdc_requests_files.py](scripts/gdc_requests_files.py)
  * [gdc_requests_legacy.py](scripts/gdc_requests_legacy.py)
  * [gdc_write_anno.py](scripts/gdc_write_anno.py)
3. generate a manifest for bulk download with `gdc-client` and a pheno file for parsing

## Parse and organize available files

To organize the files available, samples were grouped by subject into tumor-normal pairs.

### Original WRP solution
1. generate file with following fields: sample_id | sample_type | file_name | file_id
2. list normal/tumor pairs (in that order, i.e. normal first)
  * unclear how to manage duplicated samples (i.e. multiple tumors or normals for single subject)
3. executed with the following code:

```
parse_tcga_info.py hall_tcga_t10b10.tab > hall_tcga_t10b10.file_info
```

### CHD solution
1. modify parse_tcga_info.py to take as input pheno.tsv and generate samples.file_info
2. for duplicates, take sample with max sequences (i.e. reads)
  * [gdc_parse_info.py](scripts/gdc_parse_info.py)
3. output is `coad-read.file_info`

## Download ...

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
