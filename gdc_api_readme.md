# GDC API
* [Data model](https://docs.gdc.cancer.gov/Data/Data_Model/GDC_Data_Model/)
* [API guide](https://docs.gdc.cancer.gov/API/Users_Guide/Getting_Started/)

## endpoints
```
https://api.gdc.cancer.gov/<endpoint>
https://api.gdc.cancer.gov/legacy/<endpoint>

https://api.gdc.cancer.gov/annotations
https://api.gdc.cancer.gov/files
https://api.gdc.cancer.gov/projects
```

## properties
```
program.name
project.code
project_id = 'program.name-project.code'
id = uuid for entity
submitter_id = any string, unique within project, same as submitted_subject_id of study participant in dbGaP record when entity is type case
```

## types
* `project, case, demographic, sample, read_group`

## uuid
* `files, cases, samples` are examples of entities (objects) and have UUIDs

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
# get mappings
endpt = 'https://api.gdc.cancer.gov/projects/'
map = '_mapping'
response = requests.get(endpt + map)
# show in stout
print(json.dumps(response.json(), indent=2))
# write to file
with open('/scratch/chd5n/aneuploidy/mapping.txt', 'w+') as f:
    f.write(json.dumps(response.json(), indent=2))
```

## basic uuid request with curl
```
curl https://api.gdc.cancer.gov/files/0a3db3e4-8c4f-4892-b808-98c66a26c692?pretty=true
```

## basic uuid request with python
```
import requests
import json

endpt = 'https://api.gdc.cancer.gov/files/'
uuid = '0a3db3e4-8c4f-4892-b808-98c66a26c692' # file_id, use with files endpoint
response = requests.get(endpt + uuid)
print(json.dumps(response.json(), indent=2))
```

## filtered request with python
```
import requests
import json

endpt = 'https://api.gdc.cancer.gov/cases/'
filt = {
  "op": "=",
  "content": {
    "field": "cases.demographic.gender",
    "value": ["male"]
  }
}
params = {'filters': json.dumps(filt), 'fields': 'case_id'}
response = requests.get(endpt, params = params)
print(json.dumps(response.json(), indent=2))
```

## exome capture kits
```
import requests
import json

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
fields = [
  'file_id',
  'file_name',
  'analysis.metadata.read_groups.target_capture_kit_catalog_number',
  'analysis.metadata.read_groups.target_capture_kit_name',
  'analysis.metadata.read_groups.target_capture_kit_target_region',
  'analysis.metadata.read_groups.target_capture_kit_vendor'
]
fields = ','.join(fields)
params = {
  "filters": json.dumps(filt),
  "format": "JSON",
  "fields": fields,
  "size": "10"
}
response = requests.get(endpt, params = params)
print(json.dumps(response.json(), indent=2))
```

## important cases fields
* `files.downstream_analyses.output_files.msi_score, files.msi_score`

## important files fields
* `files.msi_status, files.msi_score, msi_score`

## ways to view
```
print(json.dumps(response.json(), indent=2))

with open('path/to_file.txt', 'w+') as f:
    f.write(response.content.decode('utf-8'))
```
