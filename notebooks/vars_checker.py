import pandas as pd
import json
import requests,os
from tqdm import tqdm
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

# Uncomment all the functions before running the code. Place the file in the same directory as the obs.
# csv files for KPMP and HuBMAP. There are two files in total; I am still working on the other one.

Hubmap_OPMI_onto_id = ["age", "sex", "race",
                       "azimuth_label", "predicted_label","height","weight","cause_of_death"]


def json_saver(data, output_path):
    with open(output_path, 'w') as f:
        json.dump(data, f, indent=4)


def process_datasets(datasets, output_dir):
    output_file = os.path.join(output_dir, "datasets_columns.json")
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.dirname(output_dir), exist_ok=True)
    with open(output_dir, "wb") as f:
     f.write(requests.get("https://data.bioontology.org/ontologies/OPMI/download?apikey=8b5b7825-538d-40e0-9e9e-5ab9274a9aeb&download_format=csv").content)
    columns_dict = {}

    for dataset_name, file_path in datasets.items():
        try:
            df = pd.read_csv(file_path)
            columns_dict[dataset_name] = list(df.columns)
        except Exception as e:
            print(f"Error reading {dataset_name}: {e}")

    for name, columns in columns_dict.items():
        filtered_items = [
            item for item in columns if not item.endswith("_ontology_term_id")]
        columns_dict[name] = filtered_items

    common_values = [item for item in columns_dict["KPMP SC RNAseq"]
                     if item in columns_dict["KPMP SN RNAseq"]]

    for name, columns in columns_dict.items():
        columns_dict[name] = [
            item for item in columns if item not in common_values]

    columns_dict["common-KPMP"] = common_values

    with open(output_file, "w") as f:
        json.dump(columns_dict, f, indent=4)

    print(f"Column names saved to {output_file}")


datasets = {
    'KPMP SC RNAseq': "kpmp-sc-rnaseq.obs.csv",
    'KPMP SN RNAseq': "kpmp-sn-rnaseq.obs.csv",
    'HuBMAP Left Kidney': "hubmap-LK-processed.obs.csv",
    'HuBMAP Right Kidney': "hubmap-RK-processed.obs.csv",
}

# process_datasets(datasets, "checker_output")


def check_the_onoto_label_values(value, df):
    unique_values = df[value].drop_duplicates().tolist()
    processed_values = ["Unknown" if pd.isna(v) else v for v in unique_values]
    return {value: processed_values}


def onotology_checker(datasets, output_dir, hubmap_OPMI_list):
    output_file = os.path.join(output_dir, "onotology_checker.json")
    os.makedirs(output_dir, exist_ok=True)
    onto_list = {}
    other_onto_items=["clusterClass","disease_category","assay","organism"]
    total_steps = sum([len(hubmap_OPMI_list) if name.startswith("HuBMAP") else 1
                       for name in datasets.keys()])
    progress_bar = tqdm(total=total_steps,
                        desc="Processing datasets", unit="step")

    for name, path in datasets.items():
        onto_list[name] = []
        df = pd.read_csv(path)
        if name.startswith("KPMP"):
            onot_id = df.columns[df.columns.str.endswith("_ontology_term_id")]
            onoto_label = [item[:-17] for item in onot_id]
            for item in other_onto_items:
                if item in df.columns:
                  onoto_label.append(item) 
            for item in onoto_label:
                label_values = check_the_onoto_label_values(item, df)
                onto_list[name].append(label_values)
                progress_bar.update(1)
        if hubmap_OPMI_list and name.startswith("HuBMAP"):
            for item in hubmap_OPMI_list:
                if item in df.columns:
                    label_values = check_the_onoto_label_values(item, df)
                    onto_list[name].append(label_values)
                progress_bar.update(1)

    progress_bar.close()

    with open(output_file, "w") as f:
        json.dump(onto_list, f, indent=4)
# onotology_checker(datasets, "checker_output",Hubmap_OPMI_onto_id)


def extract_fundamental_values(input_json):

    with open(input_json, "r") as f:
        data = json.load(f)

    BMI_OPMI = [
        "Body Mass Index Less Than 19",
        "Body Mass Index 19 to Less Than 21",
        "Body Mass Index 21 to Less Than 23",
        "Body Mass Index Greater Than or Equal to 23",
        "BMI greater than 40 kg per m2"
    ]
    fundamental_values = []

    for dataset, items in data.items():
        for item in items:
            for key, values in item.items():
                if isinstance(values, list):
                    fundamental_values.extend(values)

    fundamental_values = list(set(fundamental_values))
    fundamental_values.extend(BMI_OPMI)
    return fundamental_values


def csv_to_json(input_csv, output_json):
    try:
        df = pd.read_csv(input_csv)
    except FileNotFoundError as e:
        print(f"Error reading CSV file: {e}")
        return
    except pd.errors.EmptyDataError:
        print("CSV file is empty.")
        return

    key_column = df.columns[1]
    result = {}
    OPMI_in_KPMP_and_HUBmAP = {}

    for _, row in df.iterrows():
        key = row[key_column]
        row_data = row.drop(labels=key_column).to_dict()
        result[key] = row_data

    try:
        fundamental_values = extract_fundamental_values(
            r"checker_output\onotology_checker.json")
    except Exception as e:
        print(f"Error extracting fundamental values: {e}")
        return

    for item in fundamental_values:
        if isinstance(item, str):
            class_id = result.get(item, result.get(
                item.lower(), {})).get("Class ID", None)
            OPMI_in_KPMP_and_HUBmAP[item] = (
                class_id.replace("http://purl.obolibrary.org/obo/",
                                 "") if class_id else "Class ID not found"
                if class_id is not None
                else "Not exist in OPMI"
            )
        else:
            OPMI_in_KPMP_and_HUBmAP[item] = "Unknown type"

    try:
        with open(output_json, "w") as f:
            json.dump(OPMI_in_KPMP_and_HUBmAP, f, indent=4)
        print(f"JSON file successfully written to {output_json}")
    except IOError as e:
        print(f"Error writing JSON file: {e}")


csv_to_json("checker_output\OPMI.csv","checker_output\OPMI_for_KPMP_and_HUBMAP.json")


def extract_obs_columns(name, csv_file):
    df = pd.read_csv(csv_file)
    return {name: df.columns.tolist()}


def read_columns_of_file():
    output_json = {}
    for name, path in datasets.items():
        output_json.update(extract_obs_columns(name, path))
    json_saver(output_json, "checker_output\columns_for_unchange.json")


# read_columns_of_file()
