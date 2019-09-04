# GDC API
* [Data model](https://docs.gdc.cancer.gov/Data/Data_Model/GDC_Data_Model/)
* [API guide](https://docs.gdc.cancer.gov/API/Users_Guide/Getting_Started/)

## summary of coad and read using curl
```
curl 'https://api.gdc.cancer.gov/projects/TCGA-COAD?expand=summary,summary.experimental_strategies,summary.data_categories&pretty=true'
curl 'https://api.gdc.cancer.gov/projects/TCGA-READ?expand=summary,summary.experimental_strategies,summary.data_categories&pretty=true'
```

## summary of coad and read using python
```
endpt = 'https://api.gdc.cancer.gov/projects/'
project_id = 'TCGA-COAD'
params = {
  'expand': 'summary,summary.experimental_strategies,summary.data_categories'  # this cannot have spaces
}
response = requests.get(endpt + project_id, params = params)
print(json.dumps(response.json(), indent=2))

project_id = 'TCGA-READ'
response = requests.get(endpt + project_id, params = params)
print(json.dumps(response.json(), indent=2))
```

## use of mappings endpoint
```
endpt = 'https://api.gdc.cancer.gov/projects/'
project_id = '_mapping'
response = requests.get(endpt + project_id)
print(json.dumps(response.json(), indent=2))
```

## exome capture kits
```
endpt = 'https://api.gdc.cancer.gov/files/'
filt = {
  "op": "in",
  "content": {
    "field": "files.file_id",
    "value": [
      "00b7ee04-b97b-49ee-b27b-13f6fb981fec",
      "00ec5e0d-ffb4-4975-9bec-75d88a5a5411"
    ]
  }
}
fields = 'file_id,file_name,analysis.metadata.read_groups.target_capture_kit_catalog_number,analysis.metadata.read_groups.target_capture_kit_name,analysis.metadata.read_groups.target_capture_kit_target_region,analysis.metadata.read_groups.target_capture_kit_vendor'
params = {
  "filters": json.dumps(filt),
  "format": "JSON",
  "fields": fields,
  "size": "10"
}
response = requests.get(endpt, params = params)
print(json.dumps(response.json(), indent=2))
```
* [Baylor kit](https://sequencing.roche.com/en/products-solutions/by-category/target-enrichment/hybridization/seqcap-ez-hgsc-vcrome.html)
