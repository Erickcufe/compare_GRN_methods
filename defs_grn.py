import scanpy as sc
import pandas as pd
from arboreto.algo import genie3, grnboost2
from pySCENIC.utils import modules_from_adjacencies
from pySCENIC.prune import prune2df, df2regulons
from pySCENIC.aucell import aucell
from igraph import Graph
import numpy as np
import matplotlib.pyplot as plt

def infer_grn(adata):
    cell_types = adata.obs['cell_type'].unique().tolist()
    grn_results = {}

    # Load TF list from JASPAR
    tf_list = pd.read_csv('http://jaspar.genereg.net/api/v1/matrix/?tax_group=vertebrates&format=tsv', sep='\t', header=None)
    tf_list = tf_list[0].str.split('_', expand=True)[1].unique().tolist()

    for cell_type in cell_types:
        adata_subset = adata[adata.obs['cell_type'] == cell_type, :]
        expr_matrix = adata_subset.X.T

        # Run SCENIC
        adjacencies = genie3(expr_matrix, tf_names=tf_list)
        modules = modules_from_adjacencies(adjacencies, expr_matrix)
        df = prune2df(dbs=['hgnc'], modules=modules, expression_mtx=expr_matrix)
        regulons = df2regulons(df)
        auc_mtx = aucell(expr_matrix, regulons)
        grn_results[cell_type] = {'SCENIC': auc_mtx}

        # Run GENIE3
        genie3_grn = genie3(expr_matrix, tf_names=tf_list)
        grn_results[cell_type]['GENIE3'] = genie3_grn

        # Run GRNBoost2
        grnboost2_grn = grnboost2(expr_matrix, tf_names=tf_list)
        grn_results[cell_type]['GRNBoost2'] = grnboost2_grn

    return grn_results

def process_grn_results(grn_results):
    graph_results = {}

    for cell_type, methods in grn_results.items():
        graph_results[cell_type] = {}

        for method, grn in methods.items():
            g = Graph.TupleList(grn.itertuples(index=False), directed=True)
            graph_results[cell_type][method] = {
                'degree': g.degree(),
                'degree_in': g.indegree(),
                'degree_out': g.outdegree(),
                'edge_betweenness': g.edge_betweenness()
            }

    return graph_results

def calculate_rmse(graph_results):
    rmse_results = {}

    for cell_type, methods in graph_results.items():
        method_pairs = []
        rmse_values = []

        method_names = list(methods.keys())
        num_methods = len(method_names)

        for i in range(num_methods - 1):
        for j in range(i + 1, num_methods):
            method1 = method_names[i]
            method2 = method_names[j]

            degree1 = np.array(methods[method1]['degree'])
            degree2 = np.array(methods[method2]['degree'])

            rmse = np.sqrt(np.mean((degree1 - degree2) ** 2))
            method_pairs.append(f"{method1} - {method2}")
            rmse_values.append(rmse)

    rmse_results[cell_type] = {
        'method_pairs': method_pairs,
        'rmse_values': rmse_values
    }

    return rmse_results

def plot_rmse(rmse_results):
    num_cell_types = len(rmse_results)
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(15, 10))
    for idx, (cell_type, rmse_data) in enumerate(rmse_results.items()):
    ax = axes.flatten()[idx]
    ax.bar(rmse_data['method_pairs'], rmse_data['rmse_values'], color='skyblue')
    ax.set_title(f"RMSE for cell type: {cell_type}")
    ax.set_xlabel("Method pairs")
    ax.set_ylabel("RMSE")

    plt.tight_layout()
    plt.show()
