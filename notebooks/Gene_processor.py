import anndata
import requests
import json
import pandas as pd
import gzip
import shutil
import csv
import os
import scanpy as sc
from data_process import (
    common_data_across_all_daatsets,
    redundant_fields_that_need_to_remove_from_KPMP,
    extract_unique_fields_for_KPMP,
    extract_combined_fields,
    OPMI_checker,
    normalized_age,
    ontologize_BMI_to_OPMI,
    decompress_gz_to_csv,
)

# Run vars_checker.py before running data_process.py and gene_processor.py. 
# Ensure all these files are located in the same directory as the .h5ad and .csv files.


datasets = {
    "KPMP SC RNAseq": "kpmp-sc-rnaseq.h5ad",
     'KPMP SN RNAseq': "kpmp-sn-rnaseq.h5ad",
     'HuBMAP Left Kidney': "hubmap-LK-processed.h5ad",
    'HuBMAP Right Kidney': "hubmap-RK-processed.h5ad",
}

script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)
os.makedirs("main_output", exist_ok=True)

gene_name = ["SPP1"]
combined_obs_rows = []

def download_var_as_csv():
    for name, path in datasets.items():
        print("Currently working on " + name)
        adata = sc.read_h5ad(path, backed="r")
        output_file = os.path.splitext(path)[0] + ".var.csv"
        adata.var.to_csv(output_file)

def find_the_obs_by_gene_name(file_name, file_path):
    adata = sc.read_h5ad(file_path, backed="r")

    if "feature_name" in adata.var.columns:
        gene_column = "feature_name"
    elif "hugo_symbol" in adata.var.columns:
        gene_column = "hugo_symbol"
    else:
        raise ValueError(f"Neither 'feature_name' nor 'hugo_symbol' column is present in `adata.var` for {file_name}.")

    gene_in_vars = adata.var.index[adata.var[gene_column].isin(gene_name)]
    if len(gene_in_vars) == 0:
        print(f"No matching genes found in {file_name}")
        return pd.DataFrame()  
    obs_row = adata.obs[adata.X[:, adata.var.index.isin(gene_in_vars)].sum(axis=1) > 0]
    obs_row = obs_row.copy()
    obs_row["gene_name"] = ", ".join(gene_name)
    obs_row["gene_id"] = gene_in_vars[0] 
    return obs_row

def normalize_KPMP_data(row, collection):
    normalized_row = {
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
        "K-SpecimenID": row.get("SpecimenID", row.get("specimen", "Unknown")),
        "gene_name": row["gene_name"],
        "gene_id": row["gene_id"],
    }
    unique_fields = extract_unique_fields_for_KPMP()
    for value in unique_fields:
        normalized_value = row.get(value[2:], "")
        normalized_row[value] = OPMI_checker(normalized_value)
    return normalized_row

def normalize_HUBMAP_data(row, collection):
    normalized_row = {
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
        "K-SpecimenID": "",
        "gene_name": row["gene_name"],
        "gene_id": row["gene_id"],
    }
    unique_fields = extract_unique_fields_for_KPMP()
    for value in unique_fields:
        if value[2:] == "bmi":
            normalized_row["H-bmi"] = ontologize_BMI_to_OPMI(row.get(value[2:], ""))
            continue
        normalized_value = row.get(value[2:], "")
        normalized_row[value] = OPMI_checker(normalized_value)
    return normalized_row

def combine_the_obs_with_SPP1_for_kpmp_and_hubmap(output_name_and_path):
    with gzip.open(output_name_and_path + ".csv.gz", "wt", compresslevel=9, newline="") as csvfile:
        all_fields = extract_combined_fields()
        all_fields.append("gene_id")
        all_fields.append("gene_name")
        writer = csv.DictWriter(csvfile, fieldnames=all_fields)
        writer.writeheader()
        print(writer.fieldnames)
        for collection, local_path in datasets.items():
            obs_csv = find_the_obs_by_gene_name(collection, local_path)
            print("Currently working on " + collection)
            for _, row in obs_csv.iterrows():
                if collection.startswith("KPMP"):
                    normalized_row = normalize_KPMP_data(row, collection)
                elif collection.startswith("HuBMAP"):
                    normalized_row = normalize_HUBMAP_data(row, collection)
                writer.writerow(normalized_row)

combine_the_obs_with_SPP1_for_kpmp_and_hubmap("main_output/obs_with_SPP1_gene")
decompress_gz_to_csv("main_output/obs_with_SPP1_gene.csv.gz", "main_output/obs_with_SPP1_gene.csv")
