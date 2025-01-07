import anndata
import requests
import json
import pandas as pd
import gzip
import shutil
import csv
import os

script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)
# Set the working directory to the location of the script
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)
os.makedirs("main_output", exist_ok=True)
with open("checker_output\OPMI_for_KPMP_and_HUBMAP.json", "r") as f:
    OPMI_data = json.load(f)

with open("checker_output\datasets_columns.json", "r") as f:
    fields = json.load(f)


def h5ad_obs_to_csv(input_h5ad):
    output_csv = input_h5ad.replace('.h5ad', '.obs.csv')
    x = anndata.read_h5ad(input_h5ad, backed='r')
    x.obs.to_csv(output_csv)


# datasets = {
#     'KPMP SC RNAseq': 'kpmp-sc-rnaseq.h5ad',
#     'KPMP SN RNAseq': 'kpmp-sn-rnaseq.h5ad',
#     'HuBMAP Left Kidney': 'hubmap-LK-processed.h5ad'
# }

with open("checker_output\OPMI_for_KPMP_and_HUBMAP.json", "r") as file:
    data = json.load(file)


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


def prefix_for_all_fields(KPMP_SC, KPMP_SN, HUBMAP, ALL_fields):
    processed_extra_var_id_list = []
    for i in ALL_fields:
        if i in KPMP_SC or KPMP_SN:
            processed_extra_var_id_list.append("K-"+i)
        if i in HUBMAP:
            processed_extra_var_id_list.append("H-"+i)
        else:
            continue
    return processed_extra_var_id_list

def ontologize_BMI_to_OPMI(value):
    BMI=float(value)
    if BMI<=19:
        return "NCIT_C138932"
    elif 19<BMI<=21:
        return "NCIT_C138933"
    elif 21<BMI<=23:
        return "NCIT_C138934"
    elif 23<=BMI:
        return "NCIT_C138935"


common_data_across_all_daatsets = ["consortium", "collection", "dataset_id",
                                   "as_id", "cl_id", "cl_label", "gene_count", "age", "sex", "race", "disease"]
# details can be found in https://docs.google.com/spreadsheets/d/18eyNBSIRjdXsr5Om7gXUB29ZnNQb1daBWWTAao8iQU8/edit?gid=1469815304#gid=1469815304
redundant_fields_that_need_to_remove_from_KPMP = ["Race", "Age_binned"]
# Race from SC,Age_binned for all
KPMP_columns_that_nedded_to_merge = {
    "SpecimenID": {"SC": "SpecimenID", "SN": "specimen"}}


def OPMI_checker(value):
    normalized_OPMI_data = {key.lower(): val for key, val in OPMI_data.items()}
    
    lower_value = value.lower() if isinstance(value, str) else value

    if isinstance(lower_value, str) and lower_value in normalized_OPMI_data:
        if normalized_OPMI_data[lower_value] != "Not exist in OPMI":
            return normalize_category(normalized_OPMI_data[lower_value])
        else:
            return value
    else:
        return value



def extract_unique_fields_for_KPMP():
    common_value_needed_to_remove = [
        "tissue", "cell_type", "nCount_RNA", "development_stage", "sex",
        "self_reported_ethnicity", "disease", "Race", "Age_binned", "LibraryID",
        "library_id", "uuid", "age", "predicted_CLID", "predicted_label",
        "n_genes", "race", "azimuth_id","specimen","SpecimenID","diseasetype","dataSource"
    ]

    combined_fields = (
        fields.get("KPMP SC RNAseq", []) +
        fields.get("KPMP Sn RNAseq", []) +
        fields.get("common-KPMP", []) +
        fields.get("HuBMAP Left Kidney", [])
    )

    all_fields_without_common = [x for x in combined_fields if x not in common_value_needed_to_remove]

    unique_fields = ["K-SpecimenID"]
    for i in all_fields_without_common:
        if i in fields.get("KPMP SC RNAseq", []) or i in fields.get("KPMP Sn RNAseq", []) or i in fields.get("common-KPMP", []):
            unique_fields.append(f"K-{i}")
        else:
            unique_fields.append(f"H-{i}")
    return unique_fields


def extract_combined_fields():
    return extract_unique_fields_for_KPMP() + [
        "consortium",
        "collection",
        "dataset_id",
        "as_id",
        "cl_id",
        "gene_count",
        "age",
        "sex",
        "race",
        "disease"
    ]


def normalize_KPMP_data(row, collection):
    normalized_row={
        "consortium": "KPMP",
        "collection": collection,
        "dataset_id": row.get("LibraryID", row.get("library_id", "Unknown")),
        "as_id": OPMI_checker(row["tissue"]),
        "cl_id": OPMI_checker(row["cell_type"]),
        "gene_count": int(float(row["nCount_RNA"])),
        "age": normalized_age(row["Age_binned"]),
        "sex": OPMI_checker(row["sex"].title()),
        "race": OPMI_checker(row["self_reported_ethnicity"]),
        "disease": OPMI_checker(row["disease"]),
        "K-SpecimenID":row.get("SpecimenID", row.get("specimen", "Unknown"))
    }
    unique_fields = extract_unique_fields_for_KPMP()
    for value in unique_fields:
     normalized_value = row.get(value[2:], "")
     normalized_row[value] = OPMI_checker(normalized_value)

    return normalized_row

def normalize_HUBMAP_data(row, collection):
    normalized_row={
        "consortium": "HuBMAP",
        "collection": collection,
        "dataset_id": row["uuid"],
        "as_id": "kidney",
        "cl_id": OPMI_checker(row["predicted_label"]),
        "gene_count": int(float(row["n_genes"])),
        "age": normalized_age(row["age"]),
        "sex": OPMI_checker(row["sex"]),
        "race": OPMI_checker(row["race"]),
        "disease": "normal",
        "K-SpecimenID":""
    }
    unique_fields = extract_unique_fields_for_KPMP()
    for value in unique_fields:
     if value[2:]=="bmi":
         normalized_row["H-bmi"]=ontologize_BMI_to_OPMI(row.get(value[2:], ""))
         continue
     else:
      normalized_value = row.get(value[2:], "")
      normalized_row[value] = OPMI_checker(normalized_value)

    return normalized_row

def decompress_gz_to_csv(input_gz_path, output_csv_path):
    with gzip.open(input_gz_path, 'rt') as gz_file:  
        with open(output_csv_path, 'w') as csv_file:  
            shutil.copyfileobj(gz_file, csv_file) 

datasets = {
    'KPMP SC RNAseq': "kpmp-sc-rnaseq.obs.csv",
    'KPMP SN RNAseq': "kpmp-sn-rnaseq.obs.csv",
    'HuBMAP Left Kidney': "hubmap-LK-processed.obs.csv",
}
def combine_the_kpmp_and_hubmap(output_name_and_path):
    with gzip.open(output_name_and_path+".csv.gz", 'wt', compresslevel=9, newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=extract_combined_fields())
        writer.writeheader()
        print(writer.fieldnames)
        for collection, local_path in datasets.items():
            obs_csv = pd.read_csv(local_path)
            print("currently working on"+collection)
            for _, row in obs_csv.iterrows():
                if collection.startswith("KPMP"):
                 normalized_row = normalize_KPMP_data(row, collection)
                if collection.startswith("HuBMAP"):
                 normalized_row=normalize_HUBMAP_data(row, collection)
                writer.writerow(normalized_row)
# combine_the_kpmp_and_hubmap("main_output/KPMP-normalized-obs")
# decompress_gz_to_csv("main_output/KPMP-normalized-obs.csv.gz","main_output/KPMP-normalized-obs.csv")
