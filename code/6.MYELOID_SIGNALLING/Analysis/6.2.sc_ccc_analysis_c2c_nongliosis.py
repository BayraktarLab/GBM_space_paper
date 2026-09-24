import os
import warnings
import logging

import pandas as pd
import numpy as np

# Visualization
import matplotlib.pyplot as plt
import seaborn as sns
from plotnine import (
    ggplot, geom_point, aes, theme, element_text, facet_grid
)

# Single-cell analysis
import scanpy as sc
import cell2cell as c2c
import liana as li
import decoupler as dc
from liana.mt import rank_aggregate
from liana.method import (
    singlecellsignalr, connectome, cellphonedb, natmi,
    logfc, cellchat, geometric_mean
)

# Local utilities
from plots import *
from utils import *


warnings.filterwarnings('ignore')
logging.disable(logging.WARNING)

use_gpu = True
if use_gpu:
    import tensorly as tl
    tl.set_backend('pytorch')


# Argument parser
parser = argparse.ArgumentParser(description="Load an .h5ad file with Scanpy")
parser.add_argument("filepath", type=str, help="Path to the .h5ad file")

# Parse arguments
args = parser.parse_args()

# Read the integrated GBM anndata
adata = sc.read_h5ad(args.filepath)


CCI_comparisons_TME = {
    # Tumour core
    'Immune 1':[
        'Resident-Microglia TAMs', 
        'Astrocyte-like TAMs',
        'Pro-inflammatory TAMs', 
        'Interferon TAMs',
        'CD8+ T cells (cytotoxic)', 
        'Naïve T cells'
        ],
    'Immune 2':[
        'Proliferatve TAMs', 'RTN1+ TAMs', 'Anti-inflammatory TAMs',
        'Monocytes', 'Dendritic cells', 'T reg',
        'Proliferatve T cells', 'NK cells 1', 'NK cells 2'
        ],
    'Vasculature':[
        'Pericytes',
        'Endothelial cells', 
        'VLMC',
        'Resident BAM TAMs',
        'CD4+ TEM cells'
        ],
    'Gliosis':[
        'Angiogenic TAMs', 
        'Stress-response TAMs', 
        'T reg', 
        'Naïve T cells'
    ]
}

CCI_comparisons_GBM = {
    'sp_batchAT12': {'Non-proliferative 1': ['AC progenitor-like 2',
   'NPC-neuronal-like 1',
   'NPC-neuronal-like 4',
   'OPC-like 1',
   'AC-gliosis-like 3',
   'Hypoxic 2'],
  'Proliferative': ['AC progenitor-like 1',
   'AC progenitor-like 3',
   'OPC-like 2',
   'Proliferative AC-OPC-like',
   'OPC-NPC-like 1',
   'OPC-neuronal-like',
   'OPC-NPC-like 3',
   'OPC-like 3',
   'NPC-neuronal-like 3',
   'Proliferative NPC-OPC-like',
   'Hypoxic 2',
   'NPC-neuronal-like 5',
   'NPC-neuronal-like 2',
   'AC progenitor-like 4'],
  'Gliosis': ['AC-gliosis-like 1', 'OPC-like 1', 'Gliosis-like', 'Hypoxic 1']},
 'sp_batchAT6': {'Non-proliferative 1': ['AC-gliosis-like 4'],
  'Proliferative': ['AC progenitor-like 1',
   'Proliferative AC-OPC-like',
   'OPC-like 3',
   'Proliferative NPC-OPC-like',
   'AC progenitor-like 3',
   'OPC-like 2',
   'OPC-NPC-like 3',
   'AC progenitor-like 2',
   'AC progenitor-like 4',
   'NPC-neuronal-like 2',
   'AC-gliosis-like 3',
   'OPC-neuronal-like',
   'OPC-NPC-like 2',
   'OPC-NPC-like 1',
   'NPC-neuronal-like 4',
   'NPC-neuronal-like 5',
   'NPC-neuronal-like 3',
   'OPC-like 1',
   'Proliferative nIPC-like',
   'OPC-like 4'],
  'Gliosis': ['AC-gliosis-like 1',
   'AC-gliosis-like 2',
   'Gliosis-like',
   'Hypoxic 1']},
 'sp_batchAT14': {'Non-proliferative 1': ['AC-gliosis-like 3',
   'AC progenitor-like 2'],
  'Proliferative': ['Proliferative NPC-OPC-like',
   'NPC-neuronal-like 4',
   'OPC-like 1',
   'AC progenitor-like 4',
   'OPC-NPC-like 3',
   'NPC-neuronal-like 1',
   'Proliferative nIPC-like',
   'OPC-NPC-like 1',
   'AC progenitor-like 1',
   'AC-gliosis-like 2',
   'NPC-neuronal-like 3',
   'Hypoxic 2',
   'Proliferative AC-OPC-like',
   'OPC-neuronal-like',
   'OPC-like 3',
   'NPC-neuronal-like 5',
   'OPC-NPC-like 2',
   'OPC-like 2',
   'AC progenitor-like 2',
   'NPC-neuronal-like 2',
   'AC-gliosis-like 1',
   'OPC-like 4'],
  'Gliosis': ['AC-gliosis-like 1', 'Gliosis-like', 'Hypoxic 1']},
 'sp_batchAT10': {'Non-proliferative 1': ['AC progenitor-like 2',
   'AC-gliosis-like 3',
   'AC progenitor-like 4',
   'OPC-NPC-like 1',
   'OPC-like 2',
   'OPC-NPC-like 3',
   'OPC-NPC-like 2',
   'OPC-like 1',
   'OPC-like 5',
   'OPC-like 3',
   'OPC-like 4'],
  'Proliferative': ['OPC-NPC-like 1',
   'OPC-NPC-like 2',
   'OPC-like 2',
   'OPC-NPC-like 3',
   'Proliferative NPC-OPC-like',
   'Proliferative AC-OPC-like',
   'NPC-neuronal-like 4',
   'NPC-neuronal-like 2'],
  'Gliosis': ['AC-gliosis-like 1', 'Gliosis-like', 'Hypoxic 1', 'Hypoxic 2']},
 'sp_batchAT4': {'Non-proliferative 1': ['AC progenitor-like 2',
   'AC-gliosis-like 2'],
  'Proliferative': ['OPC-NPC-like 3',
   'OPC-like 3',
   'Proliferative NPC-OPC-like',
   'OPC-like 2',
   'OPC-NPC-like 1',
   'OPC-like 5',
   'NPC-neuronal-like 1',
   'OPC-NPC-like 2',
   'OPC-like 1',
   'AC progenitor-like 1',
   'OPC-neuronal-like',
   'Proliferative AC-OPC-like',
   'AC-gliosis-like 3',
   'NPC-neuronal-like 3',
   'NPC-neuronal-like 2',
   'OPC-like 4',
   'NPC-neuronal-like 5',
   'AC progenitor-like 4',
   'AC progenitor-like 3'],
  'Gliosis': ['Gliosis-like',
   'AC-gliosis-like 1',
   'Proliferative nIPC-like',
   'Hypoxic 1']},
 'sp_batchAT5': {'Non-proliferative 1': ['AC progenitor-like 2'],
  'Proliferative': ['Proliferative NPC-OPC-like',
   'Proliferative AC-OPC-like',
   'OPC-NPC-like 1',
   'Proliferative nIPC-like',
   'OPC-NPC-like 3',
   'OPC-like 3',
   'NPC-neuronal-like 3',
   'NPC-neuronal-like 5',
   'NPC-neuronal-like 4',
   'OPC-NPC-like 2',
   'AC-gliosis-like 2',
   'OPC-like 2',
   'OPC-like 4',
   'AC-gliosis-like 3'],
  'Gliosis': ['AC-gliosis-like 1',
   'AC progenitor-like 3',
   'AC progenitor-like 4',
   'AC-gliosis-like 3',
   'AC-gliosis-like 2',
   'Gliosis-like',
   'Hypoxic 2',
   'Hypoxic 1']},
 'sp_batchAT15': 
 
 {
     'Non-proliferative 1': ['AC progenitor-like 2',
                             'AC progenitor-like 4'],
      'Proliferative': ['OPC-like 3',
                      'NPC-neuronal-like 3',
                      'OPC-NPC-like 1',
                      'OPC-NPC-like 2',
                      'AC progenitor-like 1',
                      'AC progenitor-like 3',
                      'OPC-NPC-like 3',
                      'AC-gliosis-like 3',
                      'Proliferative NPC-OPC-like',
                      'OPC-like 4',
                      'Proliferative nIPC-like',
                      'OPC-like 5',
                      'NPC-neuronal-like 5',
                      'OPC-like 1',
                      'Proliferative AC-OPC-like',
                      'OPC-like 2',
                      'AC-gliosis-like 2',
                      'AC progenitor-like 4',
                      'NPC-neuronal-like 4',
                      'AC-gliosis-like 1',
                      'Hypoxic 2'],
       'Gliosis': 
                  ['Gliosis-like', 
                   'AC-gliosis-like 1', 
                   'Hypoxic 1']}
                   
                   }

analysis_name = "c2c"
donor_id_list = ['AT4', 'AT5',  'AT6', 'AT10', 'AT12','AT14', 'AT15']
tme_ligand_niches = ["Immune 1", "Immune 2", "Vasculature"]
gbm_receptor_niches = ["Non-proliferative 1", "Proliferative"]

str_lg = "_".join(tme_ligand_niches)
str_rec = "_".join(gbm_receptor_niches)
adata = adata[adata.obs.donor_id.isin(donor_id_list),:]
adata.obs["spatial_context"] = "other"

ligand_cell_states = []
for ln in tme_ligand_niches:
    for donor_id in donor_id_list:
        for col in CCI_comparisons_TME[ln]:
            if col in adata.obs["annotation_granular_vasc_coarse"].cat.categories:
                ligand_cell_states.append(col)
                adata.obs.loc[(adata.obs["annotation_granular_vasc_coarse"]==col) & (adata.obs["donor_id"]==donor_id), "spatial_context"] = f"{donor_id}_TME"



rec_cell_states = []
for ln in gbm_receptor_niches:
    
    for donor_id in donor_id_list:
        for col in CCI_comparisons_GBM[f'sp_batch{donor_id}'][ln]:
            if col in adata.obs["annotation_granular"].cat.categories:
                adata.obs.loc[(adata.obs["annotation_granular_vasc_coarse"]==col) & (adata.obs["donor_id"]==donor_id), "spatial_context"] =  f"{donor_id}_MLG"
                rec_cell_states.append(col)


groupby_pairs = itertools.product(ligand_cell_states, rec_cell_states)
groupby_pairs = pd.DataFrame(groupby_pairs, columns=['source', 'target'])
adata = adata[~(adata.obs["spatial_context"]=="other"),:]


adata.obs["spatial_context"] = adata.obs["spatial_context"].astype("category")

adata.obs["spatial_context"]

sample_key = "donor_id"
li.mt.rank_aggregate.by_sample(adata, 
                        groupby = 'annotation_granular_vasc_coarse', 
                        sample_key= sample_key,
                        groupby_pairs=groupby_pairs, 
                        expr_prop=0.05, 
                        n_jobs=40,
                        n_perms=1000,
                        verbose=True,
                        return_all_lrs=False,
                        use_raw=False)

tensor = li.multi.to_tensor_c2c(liana_res=liana_res, # LIANA's dataframe containing results
                                sample_key=sample_key, # Column name of the samples
                                source_key='source', # Column name of the sender cells
                                target_key='target', # Column name of the receiver cells
                                ligand_key='ligand_complex', # Column name of the ligands
                                receptor_key='receptor_complex', # Column name of the receptors
                                score_key='magnitude_rank', # Column name of the communication scores to use
                                inverse_fun=lambda x: 1 - x, # Transformation function
                                how='outer', # What to include across all samples
                                outer_fraction=1/3., # Fraction of samples as threshold to include cells and LR pairs.
                                context_order=liana_res[sample_key].unique(), # Order to store the contexts in the tensor
                               )

from collections import defaultdict

element_dict = defaultdict(lambda: 'Unknown')
context_dict = element_dict.copy()

# initialise_context dict
context_dict = dict.fromkeys(donor_id_list)
# assign each to donor as we interested shared pattern at the donor level
for key, _ in context_dict.keys():
	context_dict[key] = key
dimensions_dict = [context_dict, None, None, None]
meta_tensor = c2c.tensor.generate_tensor_metadata(interaction_tensor=tensor,
                                                  metadata_dicts=dimensions_dict,
                                                  fill_with_order_elements=True
                                                 )

c2c.analysis.run_tensor_cell2cell_pipeline(tensor,
                                           meta_tensor,
                                           rank=None, # Number of factors to perform the factorization. If None, it is automatically determined by an elbow analysis
                                           tf_optimization='robust', # To define how robust we want the analysis to be.
                                           upper_rank=15,
                                           random_state=0, # Random seed for reproducibility
                                           tf_svd='numpy_svd', # Type of SVD to use if the initialization is 'svd'
                                           device='cuda', # Device to use. If using GPU and PyTorch, use 'cuda'. For CPU use 'cpu'
                                           output_folder=os.path.join("../../data/out_data/sc_ccc", analysis_name), # Whether to save the figures in files. If so, a folder pathname must be passed
                                          )

import pickle

file = open(os.path.join("../../data/out_data/sc_ccc", analysis_name, f"meta_tensor_{str_lg}-VS-{str_rec}.pickle"), 'wb')
pickle.dump(meta_tensor, file)
file.close()


file = open(os.path.join("../../data/out_data/sc_ccc", analysis_name, f"tensor_{str_lg}-VS-{str_rec}.pickle"), 'wb')
# dump information to that file
pickle.dump(tensor, file)
# close the file
file.close()

# %%
factors, axes = c2c.plotting.tensor_factors_plot(interaction_tensor=tensor,
                                                 metadata = meta_tensor, # This is the metadata for each dimension
                                                 sample_col='Element',
                                                 group_col='Category',
                                                 meta_cmaps = ['tab10', 'Dark2_r', 'tab20', 'tab20'],
                                                 fontsize=10, # Font size of the figures generated
                                                 )

factors = c2c.io.load_tensor_factors(os.path.join("../../data/out_data/sc_ccc", analysis_name,'Loadings.xlsx'))

condition_colors = c2c.plotting.aesthetics.get_colors_from_labels(donor_id_list,
                                                                  cmap='plasma')
# Map these colors to each sample name
color_dict = {k : condition_colors[v] for k, v in context_dict.items()}

# Generate a dataframe used as input for the clustermap
col_colors = pd.Series(color_dict)
col_colors = col_colors.to_frame()
col_colors.columns = ['Donor id']


sample_cm = c2c.plotting.loading_clustermap(factors['Contexts'],
                                            use_zscore=False, # Whether standardizing the loadings across factors
                                            col_colors=col_colors, # Change this to color by other properties
                                            figsize=(16, 6),
                                            dendrogram_ratio=0.3,
                                            cbar_fontsize=12,
                                            tick_fontsize=14,
                                            filename= os.path.join("../../data/out_data/sc_ccc", analysis_name, f'Clustermap-Contexts_by_{sample_key}_{str_lg}-VS-{str_rec}.pdf')
                                           )


plt.sca(sample_cm.ax_heatmap)
legend = c2c.plotting.aesthetics.generate_legend(color_dict=  condition_colors  ,
                                                 bbox_to_anchor=(1.1, 0.5), # Position of the legend (X, Y)
                                                 title='Donor'
                                                )

lr_cm = c2c.plotting.loading_clustermap(factors['Ligand-Receptor Pairs'],
                                        loading_threshold=0.08, # To consider only top LRs
                                        use_zscore=True, # Whether standardizing the loadings across factors
                                        figsize=(28, 8),
                                        filename=os.path.join("../../data/out_data/sc_ccc", analysis_name, f'Clustermap-LRs_{str_lg}_{str_rec}.pdf'),
                                        row_cluster=False # To avoid clustering of factors
                                       )



