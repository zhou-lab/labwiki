cat <<EOF >1
{
    "filters":{
        "op":"in",
        "content":{
            "field":"files.file_id",
            "value":[
EOF

csvtk csv2tab Metadata-TCGA-Kraken-17625-Samples.csv | cut -f2 | awk 'NR>1{if(NR>2) printf(",\n"); printf("\""tolower($1)"\"");}' >>1

cat <<EOF >>1
]
        }
    },
    "format":"TSV",
    "fields":"file_id,file_name,cases.submitter_id,cases.case_id,data_category,data_type,cases.samples.tumor_descriptor,cases.samples.tissue_type,cases.samples.sample_type,cases.samples.submitter_id,cases.samples.sample_id,analysis.workflow_type,cases.project.project_id,cases.samples.portions.analytes.aliquots.aliquot_id,cases.samples.portions.analytes.aliquots.submitter_id",
    "size":"100000"
}
EOF

## see https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/
curl --request POST --header "Content-Type: application/json" --data @1  'https://api.gdc.cancer.gov/legacy/files' >1out
