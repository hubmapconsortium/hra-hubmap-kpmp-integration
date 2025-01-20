import anndata
import requests
import json
import pandas as pd
import gzip
import shutil
import csv
import os
import scanpy as sc
from scipy.sparse import issparse
from difflib import SequenceMatcher
import itertools
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.metrics.pairwise import cosine_similarity

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
# there are some problem with this file, comment last two rows in "data_process.py" before runing this file
# looks like the function is globle, will fix it

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
        raise ValueError(
            f"Neither 'feature_name' nor 'hugo_symbol' column is present in `adata.var` for {file_name}.")

    gene_in_vars = adata.var.index[adata.var[gene_column].isin(gene_name)]
    if len(gene_in_vars) == 0:
        print(f"No matching genes found in {file_name}")
        return pd.DataFrame()
    obs_row = adata.obs[adata.X[:, adata.var.index.isin(
        gene_in_vars)].sum(axis=1) > 0]
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
            normalized_row["H-bmi"] = ontologize_BMI_to_OPMI(
                row.get(value[2:], ""))
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
decompress_gz_to_csv("main_output/obs_with_SPP1_gene.csv.gz",
                     "main_output/obs_with_SPP1_gene.csv")


def process_and_write_in_chunks(name, path, output_path, chunk_size=10000):
    adata = sc.read_h5ad(path, backed="r")
    obs_index = adata.obs.index
    var_index = adata.var.index
    total_cells = len(obs_index)

    collection = "Unknown"
    if name.startswith("KPMP SC"):
        collection = "KPMP SC"
    elif name.startswith("KPMP SN"):
        collection = "KPMP SN"
    elif name.startswith("HuBMAP Left Kidney"):
        collection = "HuBMAP Left Kidney"
    elif name.startswith("HuBMAP Right Kidney"):
        collection = "HuBMAP Right Kidney"

    for start in range(0, total_cells, chunk_size):
        end = min(start + chunk_size, total_cells)
        print(f"Processing cells {start} to {end} for dataset: {name}")

        chunk = adata.X[start:end]
        if issparse(chunk):
            chunk = chunk.toarray()

        df_chunk = pd.DataFrame(
            chunk, index=obs_index[start:end], columns=var_index)
        df_chunk["collection"] = collection

        df_chunk.to_csv(output_path, mode="a",
                        header=not (start > 0), index=False)


def process_datasets(datasets, output_file):
    open(output_file, "w").close()

    for name, path in datasets.items():
        print(f"Starting dataset: {name}")
        process_and_write_in_chunks(name, path, output_file, chunk_size=10000)


# output_file = "main_output/cell_column_gene_row.csv"
# process_datasets(datasets, output_file)
# too big to process (need more than 50 GB RAM
filters = ["normal", "CKD", "AKI"]


def replicated_OPMI_onto_checker(name):
    return {"CKD": "Mondo_0005300", "AKI": "Mondo_0002492"}.get(name, name)


def counting_number_of_SPP1_in_different_cell_type(
    path_of_obs_with_SPP1_gene_csv_file,
    gene_type_column_name,
    output_dir, disease_filter
):
    os.makedirs(output_dir, exist_ok=True)
    obs_data = pd.read_csv(path_of_obs_with_SPP1_gene_csv_file)
    print("Available columns in the dataset:", obs_data.columns)
    if gene_type_column_name not in obs_data.columns:
        raise KeyError(
            f"Column '{gene_type_column_name}' not found in the dataset.")
    for item in disease_filter:
        disease_code = replicated_OPMI_onto_checker(item)
        filtered_data = obs_data[obs_data["disease"] == disease_code]
        print(
            f"Processing disease: {item}, Filtered rows: {len(filtered_data)}")
        value_dict = filtered_data[gene_type_column_name].value_counts(
        ).to_dict()
        output_file = os.path.join(output_dir, f"SPP1_data_for_{item}.json")
        with open(output_file, "w") as json_file:
            json.dump(value_dict, json_file, indent=4)


counting_number_of_SPP1_in_different_cell_type(
    "main_output/obs_with_SPP1_gene.csv",
    "cl_id",
    "main_output/counting_number_of_SPP1_in_different_cell_type", filters
)


def create_cosine_similarity_matrix(file_path, column_name, output_csv_path):
    obs_data = pd.read_csv(file_path)
    unique_cell_types = obs_data[column_name].unique()
    vectorizer = CountVectorizer().fit_transform(unique_cell_types)
    vectors = vectorizer.toarray()
    cosine_sim_matrix = cosine_similarity(vectors)
    similarity_df = pd.DataFrame(
        cosine_sim_matrix,
        index=unique_cell_types,
        columns=unique_cell_types
    )

    similarity_df.to_csv(output_csv_path)
    print(f"Cosine similarity matrix saved to {output_csv_path}")


create_cosine_similarity_matrix(
    "main_output/obs_with_SPP1_gene.csv", "cl_id", "main_output/cell_type_cosine_similarity.csv")
