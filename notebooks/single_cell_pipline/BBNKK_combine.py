import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np
import pooch
from scipy.sparse import csr_matrix
from scipy.io import mmwrite
import re
import anndata
import scvi

saving_h5ad='/scratch/yongqunh_root/yongqunh0/wruopeng/raw_combined_output/pipline/data_integration/ingest_method'
sc.settings.figdir = "/scratch/yongqunh_root/yongqunh0/wruopeng/raw_combined_output/pipline/data_integration/PP_subset/"
if not os.path.exists(sc.settings.figdir):
    os.makedirs(sc.settings.figdir)
sc.settings.autosave = True

sc.logging.print_versions()
sc.set_figure_params(facecolor="white", figsize=(10, 10))
sc.settings.verbosity = 3

HuBMAP_LK_raw=sc.read('/scratch/yongqunh_root/yongqunh0/wruopeng/HuBMAP_KPMP/standard_pipline/HuBMAP_LK_Harmony_pre_proccessed.h5ad')
HuBMAP_RK_raw=sc.read('/scratch/yongqunh_root/yongqunh0/wruopeng/HuBMAP_KPMP/standard_pipline/HuBMAP_RK_Harmony_pre_proccessed.h5ad')
KPMP_SC_raw=sc.read('/scratch/yongqunh_root/yongqunh0/wruopeng/HuBMAP_KPMP/standard_pipline/KPMP_SC_Harmony_pre_proccessed.h5ad')

KPMP_shared = {
    "LibraryID": "dataset_id",
    "tissue_ontology_term_id": "as_id",
    "cell_type_ontology_term_id": "cl_id",
    "cell_type": "cl_label",
    # "nCount_RNA": "gene_count",
    "Age_binned": "age",
    "sex": "sex",
    "self_reported_ethnicity": "race",
    "disease": "disease",
}
KPMP_SC_raw.obs['consortium'] = 'KPMP'
KPMP_SC_raw.obs['collection'] = 'KPMP SC RNAseq'

KPMP_SC_raw.obs.rename(columns=KPMP_shared, inplace=True)

KPMP_var_shared = {
    'feature_name': 'gene_name'
}
KPMP_SC_raw.var.rename(columns=KPMP_var_shared, inplace=True)

Hubmap_shared = {
    "cell_id": "cell_id", 
    "predicted_CLID": "cl_id",
    "predicted_label": "cl_label",
    # "n_genes": "gene_count",
    "age": "age",
    "sex": "sex",
    "race": "race"
}
HubMap_var_shared = {
    "hugo_symbol": "gene_name"
}

def add_common_value(adata,RK):
    # Check if 'uuid' column exists and rename it if so.
    if 'uuid' in adata.obs.columns:
        adata.obs.rename(columns={"uuid": "dataset_id"}, inplace=True)
    else:
        adata.obs['dataset_id'] = 'currently not included in LK'
    # Set common values
    adata.obs['as_id'] = "UBERON:0002113"
    adata.obs['disease'] = "normal"
    adata.obs['consortium'] = 'HuBMAP'
    if RK == True:
        adata.obs['collection'] = 'HuBMAP Right Kidney'
    else:
        adata.obs['collection'] = 'HuBMAP Left Kidney'
    adata.obs['tissue'] = 'kidney'
    # Rename HuBMAP columns using inplace renaming
    adata.obs.rename(columns=Hubmap_shared, inplace=True)
    adata.var.rename(columns=HubMap_var_shared, inplace=True)
    adata.var.index = adata.var.index.to_series().apply(lambda x: re.sub(r'\.\d+$', '', x))
    return adata

# apply add_common_value to your HuBMAP datasets if needed:
# HuBMAP_LK_raw = add_common_value(HuBMAP_LK_raw,False)
# HuBMAP_RK_raw = add_common_value(HuBMAP_RK_raw,True)

def  remove_extra_var(adata):
    columns_to_drop = [col for col in adata.var.columns if col != 'gene_name']
    adata.var.drop(columns=columns_to_drop, inplace=True)
    return adata

# KPMP_SC_raw=remove_extra_var(KPMP_SC_raw)

def drop_unqiue_obs_columns(adata):
    keep_cols = [
    'consortium',
    'collection',
    'dataset_id',
    'as_id',
    'cl_id',
    'cl_label',
    # 'gene_count',
    'age',
    'sex',
    'race',
    'disease'
    ]
    missing_cols = [col for col in keep_cols if col not in adata.obs.columns]
    if missing_cols:
        print("Warning: The following columns are missing in adata.obs:", missing_cols)
    available_cols = [col for col in keep_cols if col in adata.obs.columns]
    adata.obs = adata.obs[available_cols]
    return adata

# HuBMAP_LK_raw = drop_unqiue_obs_columns(HuBMAP_LK_raw)
# HuBMAP_RK_raw = drop_unqiue_obs_columns(HuBMAP_RK_raw)
# KPMP_SC_raw=drop_unqiue_obs_columns(KPMP_SC_raw)


def filter_common_cl_labels(adatas):
    """
    Given a list of AnnData objects, determine the common cell type labels (cl_label)
    shared by all datasets and filter each dataset to retain only cells with those labels.
    Returns the filtered AnnData objects and prints the shared labels.
    """
    for ad in adatas:
        print(ad.obs.columns)
        print(ad.var.columns)
        print(ad.var['gene_name'].str.startswith("MT-"))
    # Get the set of unique cl_labels from the first dataset
    common_labels = set(adatas[0].obs['cl_label'].unique())
    # Intersect with cl_labels from remaining datasets
    for ad in adatas[1:]:
        common_labels &= set(ad.obs['cl_label'].unique())
    print("Common cl_label values:", common_labels)
    
    # Filter each dataset to retain only cells with a cl_label in the common set
    filtered_adatas = []
    for ad in adatas:
        ad_filtered = ad[ad.obs['cl_label'].isin(common_labels)].copy()
        # Optionally, you can reset the cl_label values to themselves
        # (this step is not necessary if they already match the common set)
        ad_filtered.obs['cl_label'] = ad_filtered.obs['cl_label'].astype('category')
        filtered_adatas.append(ad_filtered)
    return filtered_adatas

# HuBMAP_LK_raw, HuBMAP_RK_raw, KPMP_SC_raw = filter_common_cl_labels([HuBMAP_LK_raw, HuBMAP_RK_raw, KPMP_SC_raw])

def map_the_gene_name(adata):
    common_genes = KPMP_SC_raw.var.index.intersection(adata.var.index)
    gene_names = KPMP_SC_raw.var.loc[common_genes, 'gene_name'].astype(str)
    adata.var['gene_name'] = adata.var['gene_name'].astype(str)
    adata.var.loc[common_genes, 'gene_name'] = gene_names
    
    return adata

# HuBMAP_LK_raw = map_the_gene_name(HuBMAP_LK_raw)
# HuBMAP_RK_raw = map_the_gene_name(HuBMAP_RK_raw)


def filter_to_common_genes(adata_list):
    # Identify common genes across all AnnData objects
    common_genes = set(adata_list[0].var_names)
    for ad in adata_list[1:]:
        common_genes = common_genes.intersection(set(ad.var_names))
    common_genes = list(common_genes)
    
    # Subset each AnnData to only include the common genes and return a new list of AnnData objects
    filtered_adatas = [ad[:, common_genes].copy() for ad in adata_list]
    # print(f"Number of common genes kept: {len(common_genes)}")   
    return filtered_adatas

# HuBMAP_LK_raw, HuBMAP_RK_raw, KPMP_SC_raw = filter_to_common_genes([HuBMAP_LK_raw, HuBMAP_RK_raw, KPMP_SC_raw])


def QC_proccess(adata,batch,QC_params,scVI):
    adata.write_h5ad(f'/scratch/yongqunh_root/yongqunh0/wruopeng/raw_combined_output/pipline/data_integration/PP_subset/{batch}_before_proccessed.h5ad')
    adata.var_names_make_unique()
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    # mitochondrial genes, "MT-" for human, "Mt-" for mouse
    adata.var["mt"] = adata.var['gene_name'].str.startswith("MT-")
    # ribosomal genes
    # adata.var["ribo"] = adata.var['gene_name'].str.startswith(("RPS", "RPL"))
    # # hemoglobin genes
    # adata.var["hb"] = adata.var['gene_name'].str.contains("^HB[^(P)]")
    sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, .98)
    sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,save=f'_{batch}_data_distribution.png'
    )
    sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt",save=f'_{batch}_slope_distribution.png')
    adata = adata[adata.obs.n_genes_by_counts < upper_lim, :]
    adata = adata[adata.obs.pct_counts_mt < QC_params['pct_counts_mt'], :]
    # adata = adata[adata.obs['total_counts'] < QC_params['total_counts'], :]
    sc.pp.scrublet(adata)
    adata = adata[~adata.obs["predicted_doublet"], :] 
    if scVI == False:
        adata.layers["counts"] = adata.X.copy()
        # Normalizing to median total counts
        sc.pp.normalize_total(adata,target_sum=1e4)
        # Logarithmize the data
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        adata.raw = adata
        adata = adata[:, adata.var.highly_variable]
        sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
        sc.pp.scale(adata, max_value=10)

        sc.tl.pca(adata, svd_solver='arpack')
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)
        sc.tl.leiden(adata, resolution = 0.25)
        sc.tl.umap(adata)
        sc.pl.umap(adata, color="cl_label", save=f'{batch}_adter_QC_umap.png')
        adata.write_h5ad(f'/scratch/yongqunh_root/yongqunh0/wruopeng/raw_combined_output/pipline/data_integration/PP_subset/{batch}_pre_proccessed.h5ad')
        return adata
    if scVI == True:
        print('no pca pre scVI')
        adata.write_h5ad(f'{saving_h5ad}/{batch}_pre_scVI_proccessed.h5ad')
        return adata
    else:
        return adata

# HuBMAP_LK_raw=QC_proccess(HuBMAP_LK_raw,'HuBMAP_LK',{
#     'n_genes_by_counts':2000,
#     'total_counts':90000,
#     'pct_counts_mt':20
# },scVI=False)
# HuBMAP_RK_raw=QC_proccess(HuBMAP_RK_raw,'HuBMAP_RK',{
#     'n_genes_by_counts':2000,
#     'total_counts':90000,
#     'pct_counts_mt':20
# },scVI=False)
# KPMP_SC_raw=QC_proccess(KPMP_SC_raw,'KPMP_SC',{
#     'n_genes_by_counts':5000,
#     'total_counts':20000,
#     'pct_counts_mt':100
# },scVI=False)



def data_integration_BBNK(adata_list, batch_keys):
    # 3. Concatenate all datasets (note: ingest step removed)
    adata_concat = anndata.concat(adata_list, label="batch", keys=batch_keys,
    axis=0, join="inner", merge="first")
    sc.external.pp.bbknn(adata_concat, batch_key="batch")
    sc.tl.leiden(adata_concat,resolution='0.6')
    sc.tl.umap(adata_concat)
    sc.pl.umap(adata_concat, color=['batch','cl_label','disease'],save='bbnkk_LK_KPMP_umap.png')
    
    return adata_concat

def data_integration_ingest(adata_list, batch_keys, output_dir="."):
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Find common genes across all datasets
    common_genes = set(adata_list[0].var_names)
    for ad in adata_list[1:]:
        common_genes = common_genes.intersection(set(ad.var_names))
    common_genes = list(common_genes)
    
    # Subset each AnnData to only include the common genes
    for i in range(len(adata_list)):
        adata_list[i] = adata_list[i][:, common_genes]
    
    # Ingest the reference dataset's leiden clusters into all query datasets
    for ad in adata_list[1:]:
        sc.tl.ingest(ad, adata_list[0], obs='leiden')
    
    # 2. Designate the first dataset as the reference and process it
    adata_ref = adata_list[0]
    # sc.pp.pca(adata_ref)
    # sc.pp.neighbors(adata_ref)
    # sc.tl.umap(adata_ref)
    # sc.pl.umap(adata_ref, color="cl_label", save='ingest_adata_ref_umap.png', show=False)
    
    # 3. Ingest each query dataset into the reference (map cluster labels)
    for i in range(1, len(adata_list)):
        sc.tl.ingest(adata_list[i], adata_ref, obs='cl_label')
        adata_list[i].uns['cl_label_colors'] = adata_ref.uns['cl_label_colors']
        sc.pl.umap(adata_list[i], color=['cl_label'], wspace=0.5,
                     save=f'ingest__cl_label_bulk_labels_{i}.png', show=False)
    
    # 4. Concatenate all datasets, assigning the provided batch keys
    adata_concat = anndata.concat(adata_list, label="batch", keys=batch_keys,
                                  axis=0, join="inner", merge="first")
    concat_file = os.path.join(output_dir, 'ingest_adata_concat_ingest.h5ad')
    adata_concat.write_h5ad(concat_file)
    
    # Plot integrated UMAP using the concatenated AnnData object
    sc.pl.umap(adata_concat, color='leiden', frameon=False, title='Integrated', legend_loc=None)
    sc.pl.umap(adata_concat, color=['batch','cl_label'])
    
    return adata_concat

def data_integration_SCVI(adata_list, batch_keys,save_path):
    if not os.path.exists(save_path):
      os.makedirs(save_path)
    common_genes = set(adata_list[0].var_names)
    for ad in adata_list[1:]:
        common_genes = common_genes.intersection(set(ad.var_names))
    common_genes = list(common_genes)
    
    # Subset each AnnData to only include the common genes
    for i in range(len(adata_list)):
        adata_list[i] = adata_list[i][:, common_genes]

    adata = anndata.concat(adata_list, label="Sample", keys=batch_keys,axis=0, join="inner", merge="first")
    # adata.write_h5ad(f'{save_path}/Pre_SCVI.h5ad')
    sc.pp.filter_cells(adata, min_genes = 200)
    sc.pp.filter_genes(adata, min_cells = 3)
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars = ['mt'], percent_top = None, log1p = False, inplace = True)
    adata = adata[adata.obs.pct_counts_mt < 50]
    adata.layers['counts'] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum = 1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, subset = True, layer = 'counts',flavor = "seurat_v3", batch_key="Sample")
    scvi.model.SCVI.setup_anndata(adata, layer = "counts",categorical_covariate_keys=["Sample"],continuous_covariate_keys=['pct_counts_mt', 'total_counts'])
    model = scvi.model.SCVI(adata)
    model.train()
    latent = model.get_latent_representation()
    adata.obsm['X_scVI'] = latent
    adata.layers['scvi_normalized'] = model.get_normalized_expression(library_size = 1e4)
    sc.pp.neighbors(adata, use_rep = 'X_scVI')
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution = 0.5)
    sc.pl.umap(adata, color = ['leiden', 'Sample','cl_label'], frameon = False,save='_SCVI_method_umap.png')
    model.save(f"{save_path}/KPMP_HuBMAP_model", overwrite=True)
    adata.obs['age'] = adata.obs['age'].astype(str)
    return adata

def data_integration_harmony(adata_list, batch_keys):
    adata = anndata.concat(adata_list, label="Sample", keys=batch_keys,axis=0, join="inner", merge="first")
    sc.tl.pca(adata, n_comps=50, svd_solver='arpack')
    sc.external.pp.harmony_integrate(adata, key='batch', basis='X_pca')
    sc.pp.neighbors(adata, use_rep='X_pca_harmony')
    sc.tl.umap(adata)
    sc.tl.leiden(adata)
    sc.pl.umap(adata, color=['batch','cl_label','disease'],save='harmony_umap.png')
    return adata


# Combined_data_SCVI = data_integration_SCVI([KPMP_SC_raw, HuBMAP_LK_raw, HuBMAP_RK_raw],
#                            ["KPMP_SC", "HuBMAP_LK", "HuBMAP_RK"],"/scratch/yongqunh_root/yongqunh0/wruopeng/raw_combined_output/pipline/data_integration/SCVI_method")

# Combined_data_SCVI.write_h5ad('/scratch/yongqunh_root/yongqunh0/wruopeng/raw_combined_output/pipline/data_integration/SCVI_method/Combined_data_SCVI.h5ad')

# Combined_data_ingest=data_integration_ingest([
#     # KPMP_SC_raw,
#  HuBMAP_LK_raw, HuBMAP_RK_raw],
#                            [
#                             # "KPMP_SC",
#                             "HuBMAP_LK", "HuBMAP_RK"],'/scratch/yongqunh_root/yongqunh0/wruopeng/raw_combined_output/pipline/data_integration/ingest_method')
# Combined_data_ingest.write_h5ad('/scratch/yongqunh_root/yongqunh0/wruopeng/raw_combined_output/pipline/data_integration/ingest_method/Combined_data_HuBMAP_ingest.h5ad')

# Now, call the function with your three datasets:
Combined_data_BBNK = data_integration_BBNK([KPMP_SC_raw, HuBMAP_LK_raw, HuBMAP_RK_raw],
                           ["KPMP_SC", "HuBMAP_LK", "HuBMAP_RK"])

Combined_data_BBNK.write_h5ad('/scratch/yongqunh_root/wruopeng/HuBMAP_KPMP/standard_pipline/Combined_data_BBNK.h5ad')