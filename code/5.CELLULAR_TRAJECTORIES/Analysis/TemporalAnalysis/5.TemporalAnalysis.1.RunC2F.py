import scvelo as scv
import scanpy as sc
import cell2fate as c2f
import pickle as pickle
import pandas as pd
import numpy as np
from os.path import exists
import matplotlib.pyplot as plt
import torch
import random
import scvi
import argparse


data_dir = '../data/GBM/'
clone_dir = '../data/GBM/clones/'
save_dir = '../data/cell2fate/'

def main(clone_name: str ):
    # clone_name = patient + 'clone' + str(unique_clones[j])
    adata = sc.read_h5ad(clone_dir + clone_name + '_2.h5ad')
    adata = adata[adata.obs["TME_GBM_granular"].str.contains("like")]
    adata.layers['counts'] = adata.X
    adata = c2f.utils.get_training_data(adata, cells_per_cluster = 10**6, cluster_column = "TME_GBM_granular",
                                    remove_clusters = [], min_shared_counts = 20, n_var_genes= 3000)
    c2f.Cell2fate_DynamicalModel.setup_anndata(adata, spliced_label='spliced', unspliced_label='unspliced',
                                                batch_key = 'sample')
    print('adata.shape', adata.shape)
    n_modules = c2f.utils.get_max_modules(adata)
    print('n_modules', n_modules)
    from cell2fate._pyro_mixin import MyAutoHierarchicalNormalMessenger
    mod = c2f.Cell2fate_DynamicalModel(adata, n_modules = n_modules,
                                        guide_class=MyAutoHierarchicalNormalMessenger)
    patience = 50
    mod.train(batch_size = 1000, max_epochs = 5000,
                early_stopping = True, early_stopping_min_delta = 10**(-4),
                early_stopping_monitor = 'elbo_train', early_stopping_patience = patience)
    adata=mod.export_posterior_quantiles(adata, batch_size=1000)
    post_sample_means = mod.samples['post_sample_means']
    post_sample_means['var_names'] = adata.var_names
    post_sample_means['obs_names'] = adata.obs_names
    post_sample_means['ELBO'] = np.mean(np.array(mod.history['elbo_train'][-patience:]))
    post_sample_means['TME_GBM_granular'] = adata.obs['TME_GBM_granular']
    with open(save_dir + clone_name + '_2_new2.pickle', 'wb') as handle:
        pickle.dump(post_sample_means, handle, protocol=pickle.HIGHEST_PROTOCOL)
    adata.write_h5ad(save_dir + clone_name + '_2_new2.h5ad')


if __name__ == "__main__":
    parser = argparse.ArgumentParser() 
    parser.add_argument("clone_name", type=str) 
    args = parser.parse_args()
    main(args.clone_name)