import requests
import json
import pandas as pd
import gzip
import shutil


def processed_fields_list():
    updated_fields = ["consortium", "collection", "dataset_id", "cell_id", "as_id",
                      "cl_id", "cl_label", "gene_count", "age", "sex", "race", "disease"]
    hubmap_id_list = [
        "uuid",
        "hubmap_id",
        "age",
        "sex",
        "height",
        "weight",
        "bmi",
        "cause_of_death",
        "race",
        "barcode",
        "dataset",
        "azimuth_label",
        "azimuth_id",
        "predicted_CLID",
        "predicted_label",
        "cl_match_type",
        "prediction_score",
        "n_genes",
        "n_counts",
        "leiden"
    ]
    extra_id_list = ["subclass.l1",
                     "assay_ontology_term_id", "disease_ontology_term_id",
                     "suspension_type",
                     "disease_category",
                     "development_stage",
                     "organism",
                     "organism_ontology_term_id",
                     "donor_id",
                     "sex_ontology_term_id",
                     "nFeature_RNA",
                     "diabetes_history",
                     "tissue_type", "is_primary_data",
                     "tissue",
                     "hypertension",
                     "assay",
                     "development_stage_ontology_term_id",
                     "self_reported_ethnicity_ontology_term_id",
                     "observation_joinid",
                     "eGFR",
                     "percent.mt",
                     "sampletype",
                     "Race",
                     "SpecimenID",
                     "clusterClass",
                     "clusterNumber",
                     "SampleID",
                     "author_cell_type",
                     "Run",
                     "dataSource",
                     "diseasetype",
                     "experiment_id",
                     "percent.medulla",
                     "percent.er",
                     "percent.cortex",
                     "specimen",
                     "subclass.l2",
                     "region",
                     "class",
                     "hubmap_id",
                     "height",
                     "weight",
                     "bmi",
                     "cause_of_death",
                     "barcode",
                     "dataset",
                     "azimuth_label",
                     "azimuth_id",
                     "cl_match_type",
                     "prediction_score",
                     "n_counts",
                     "leiden"
                     ]
    processed_extra_id_list=[]
    for id in extra_id_list:
        if id in hubmap_id_list:
            processed_extra_id_list.append("H-"+id)
        else:
            processed_extra_id_list.append("K-"+id)
    # extra_id_list = ["N-" + var for var in extra_id_list]
    updated_fields.extend(processed_extra_id_list)
    extra_id_list=processed_extra_id_list
    return updated_fields, extra_id_list

updated_fields,extra_id_list=processed_fields_list()

def normalize_category(value):
    if not isinstance(value, str):
        value = str(value)  
    value = value.strip().lower()
    if value == "" or value == "unknown":
        return "Unknown"
    return value.title()

def normalized_age(age):
    cases = {
        "first": 0,
        "second": 10,
        "third": 20,
        "fourth": 30,
        "fifth": 40,
        "sixth": 50,
        "seventh": 60,
        "eighth": 70,
        "nineth": 80,
        "tenth": 90
    }
    # Ensure `age` is a string before splitting
    if isinstance(age, str):
        age_key = age.split(" ")[0]
        age = str(cases.get(age_key, age))
    else:
        age = str(age)  # Convert float or int to string
    
    if age != "" and age[0].isdigit():
        return f"{age[0]}0-{age[0]}9"
    else:
        return "Unknown"

def normalize_race(race):
    if race == "African American":
        return "Black or African American"
    else:
        return normalize_category(race)

_donor_cache = {}
def add_hubmap_donor_info(row):
    if row["sex"] == "Unknown" or row["age"] == "Unknown" or row["race"] == "Unknown":
        uuid = row["dataset_id"]
        if uuid in _donor_cache:
            donor = _donor_cache[uuid]
        else:
            url = f"https://entity.api.hubmapconsortium.org/datasets/{uuid}/donors"
            donor_data = json.loads(requests.get(url).content)[0]
            metadata = donor_data["metadata"]["living_donor_data"]
            donor = {}
            for md in metadata:
                key = md["grouping_concept_preferred_term"].lower()
                if key == "sex" or key == "race":
                    donor[key] = md["preferred_term"]
                elif key == "age":
                    donor[key] = md["data_value"]
                    if md["units"] == "months":
                        donor[key] = str(float(donor[key]) / 12)
                    donor[key] = normalized_age(donor[key])
            _donor_cache[uuid] = donor

        for (key, value) in donor.items():
            row[key] = value

def normalize_row_kpmp_with_unique_vars(row, collection):
    normalized_row={
        "consortium": "KPMP",
        "collection": collection,
        "dataset_id": row.get("LibraryID", row.get("library_id", "Unknown")),
        "cell_id": "",
        "as_id": row["tissue_ontology_term_id"],
        "cl_id": row["cell_type_ontology_term_id"],
        "cl_label": row["cell_type"],
        "gene_count": int(float(row["nCount_RNA"])),
        "age": normalized_age(row["Age_binned"]),
        "sex": normalize_category(row["sex"].title()),
        "race": normalize_race(row["self_reported_ethnicity"]),
        "disease": normalize_category(row["disease"])
    }
    for vars in extra_id_list:
        if vars[2:] in row:
            normalized_row[vars]=normalize_category(row[vars[2:]])
        else:
           normalized_row[vars]=""
    return normalized_row

def normalize_row_hubmap_with_unique_vars(row, collection):
    normalized_row = {
        "consortium": "HuBMAP",
        "collection": collection,
        "dataset_id": row["uuid"],
        "cell_id": row["cell_id"],
        "as_id": "UBERON:0002113",
        "cl_id": row["predicted_CLID"],
        "cl_label": row["predicted_label"],
        "gene_count": int(row["n_genes"]),
        "age": normalized_age(row["age"]),
        "sex": normalize_category(str(row["sex"]).title()),
        "race": normalize_race(row["race"]),
        "disease": "normal"
    }
    for vars in extra_id_list:
        if vars[2:] in row:
            normalized_row[vars]=normalize_category(row[vars[2:]])
        else:
           normalized_row[vars]=""
    add_hubmap_donor_info(normalized_row)
    return normalized_row

fields = ["consortium", "collection", "dataset_id", "cell_id", "as_id", "cl_id", "cl_label", "gene_count", "age", "sex", "race", "disease"]
import csv
import gzip
datasets = {
    'KPMP SC RNAseq': "E:/Git/ProfessorHe's research/kpmp-sc-rnaseq.obs.csv",
    'KPMP SN RNAseq': "E:/Git/ProfessorHe's research/kpmp-sn-rnaseq.obs.csv",
    'HuBMAP Left Kidney': "E:/Git/ProfessorHe's research/hubmap-LK-processed.obs.csv",
    'HuBMAP Right Kidney': "E:/Git/ProfessorHe's research/hubmap-RK-processed.obs.csv",
}

def combine_the_kpmp_and_hubmap(output_name_and_path):
    with gzip.open(output_name_and_path+".csv.gz", 'wt', compresslevel=9, newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=updated_fields)
        writer.writeheader()
        print(writer.fieldnames)
        for collection, local_path in datasets.items():
            obs_csv = pd.read_csv(local_path)
            print("currently working on"+collection)
            for _, row in obs_csv.iterrows():
                 if collection.startswith("HuBMAP"):
                        normalized_row = normalize_row_hubmap_with_unique_vars(row, collection)
                 else:
                        normalized_row = normalize_row_kpmp_with_unique_vars(row, collection)
                 writer.writerow(normalized_row)
combine_the_kpmp_and_hubmap("E:/Git/ProfessorHe's research/V2_combine_all_gz")