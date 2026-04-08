# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 11:02:23 2026

@author: p03565
"""

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.data import Data, DataLoader
from torch_geometric.nn import GATConv
from torch_geometric.utils import negative_sampling
from typing import List, Tuple
from tqdm import tqdm

#Import files next#
csv_file_path_1 = 'C:\\Fish Data - Dr. Camara\\Moselle\\GTN_results_First\\NB LCMM\\Hmsc + GNN+\\Extracted Files\\abundance_baseline.csv'
csv_file_path_2 = 'C:\\Fish Data - Dr. Camara\\Moselle\\GTN_results_First\\NB LCMM\\Hmsc + GNN+\\Extracted Files\\omega_baseline.csv'
csv_file_path_3 = 'C:\\Fish Data - Dr. Camara\\Moselle\\GTN_results_First\\NB LCMM\\Hmsc + GNN+\\Extracted Files\\metadata_baseline.csv'
csv_file_path_4 = 'C:\\Fish Data - Dr. Camara\\Moselle\\GTN_results_First\\NB LCMM\\Hmsc + GNN+\\Extracted Files\\traits_baseline.csv'
csv_file_path_5 = 'C:\\Fish Data - Dr. Camara\\Moselle\\GTN_results_First\\NB LCMM\\Hmsc + GNN+\\Extracted Files\\Beta Coeffs HMSC.csv'
csv_file_path_6 = 'C:\\Fish Data - Dr. Camara\\Moselle\\GTN_results_First\\NB LCMM\\Hmsc + GNN+\\Extracted Files\\species_family.xlsx'


#Load all data sets#
abundance_df = pd.read_csv(csv_file_path_1)
omega_df = pd.read_csv(csv_file_path_2)
metadata_df = pd.read_csv(csv_file_path_3)
traits_df = pd.read_csv(csv_file_path_4)
beta_df = pd.read_csv(csv_file_path_5)
species_family_df = pd.read_excel(csv_file_path_6)
#check columns 
abundance_df.columns
traits_df.columns
omega_df.columns
species_family_df.columns

omega_df
traits_df

# ==============================================
# IMPORTS & NUMPY FIX
# ==============================================
import numpy as np
import pandas as pd
import torch
from torch_geometric.data import Data
import networkx as nx
import matplotlib.pyplot as plt
from tqdm import tqdm
import imageio
import os

# Fix NumPy 2.0
if not hasattr(np, 'float_'):
    np.float_ = np.float64

# ==============================================
# PREPROCESSING
# ==============================================
def preprocess_data(abundance_df, metadata_df, species_cols):

    # Strip strings
    abundance_df['sample_id'] = abundance_df['sample_id'].astype(str).str.strip()
    metadata_df['spatiotemporal'] = metadata_df['spatiotemporal'].astype(str).str.strip()

    # Merge by row order (ensure 1:1)
    merged_df = abundance_df.copy()
    merged_df['spatiotemporal'] = metadata_df['spatiotemporal']

    # Aggregate by station_year
    agg_df = merged_df.groupby('spatiotemporal')[species_cols].mean().reset_index()
    agg_df['station_year'] = agg_df['spatiotemporal']
    agg_df['station'] = agg_df['station_year'].apply(lambda x: x.split('_')[0])
    agg_df['year'] = agg_df['station_year'].apply(lambda x: x.split('_')[1])

    print(f"✓ Rows after aggregation: {len(agg_df)}")
    print(f"✓ Unique station_year: {len(agg_df['station_year'].unique())}")

    return agg_df

# ==============================================
# GRAPH BUILDER
# ==============================================
class SectorGraphBuilder:
    def __init__(self, omega_df, traits_df, abundance_df, species_family_df, threshold=0.1):

        self.species_names = omega_df['species_code'].tolist()
        self.omega = omega_df.set_index('species_code').loc[self.species_names,self.species_names].values.astype(np.float64)
        self.traits = traits_df.set_index('species_code').loc[self.species_names].values.astype(np.float64)
        self.abundance_df = abundance_df.reset_index(drop=True)
        self.threshold = threshold

        # Family colors
        self.family_map = species_family_df.set_index('SP_Code')['FAMILY'].to_dict()
        self.families = [self.family_map.get(sp, 'Unknown') for sp in self.species_names]
        unique_families = list(set(self.families))
        cmap = plt.cm.get_cmap('tab20', len(unique_families))
        self.family_color_map = {fam: cmap(i) for i,fam in enumerate(unique_families)}

        print(f"✓ {len(self.species_names)} species loaded")

    def build_edge_index_and_weights(self):
        edges, weights = [], []
        for i in range(len(self.species_names)):
            for j in range(i+1, len(self.species_names)):
                if abs(self.omega[i,j]) > self.threshold:
                    edges += [[i,j],[j,i]]
                    weights += [self.omega[i,j]]*2
        edge_index = torch.tensor(edges, dtype=torch.long).t().contiguous()
        edge_weights = torch.tensor(weights, dtype=torch.float32).unsqueeze(1)
        return edge_index, edge_weights

    def build_graph(self, row):
        abundances = row[self.species_names].values.astype(np.float64).reshape(-1,1)
        node_features = np.hstack([self.traits, abundances])
        x = torch.tensor(node_features, dtype=torch.float32)
        edge_index, edge_weights = self.build_edge_index_and_weights()
        return Data(x=x, edge_index=edge_index, edge_attr=edge_weights)

# ==============================================
# NETWORK PLOT
# ==============================================
def plot_network(builder, row, save_path_svg, save_path_png=None):

    G = nx.Graph()
    n = len(builder.species_names)
    G.add_nodes_from(range(n))

    edge_index, edge_weights = builder.build_edge_index_and_weights()
    for i in range(edge_index.shape[1]):
        u = edge_index[0,i].item()
        v = edge_index[1,i].item()
        if u<v:
            G.add_edge(u, v, weight=edge_weights[i].item())

    abundances = row[builder.species_names].values.astype(float)
    node_sizes = 300 + 3000 * (abundances / (abundances.max() + 1e-6))
    node_colors = [builder.family_color_map[fam] for fam in builder.families]

    pos = nx.spring_layout(G, seed=42, k=2)
    fig, ax = plt.subplots(figsize=(12,10))

    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors,
                           edgecolors='black', linewidths=1, ax=ax)
    weights = [G[u][v]['weight'] for u,v in G.edges()]
    nx.draw_networkx_edges(G, pos, edge_color=weights, edge_cmap=plt.cm.RdBu,
                           edge_vmin=-1, edge_vmax=1, width=2, alpha=0.7, ax=ax)

    # labels
    for i,(x,y) in pos.items():
        ax.text(x, y, builder.species_names[i], fontsize=8, ha='center', va='center',
                bbox=dict(facecolor='white', alpha=0.6, edgecolor='none'))

    ax.set_title(f"{row['station_year']}\nNode size = abundance | Color = family | Edge = Ω",
                 fontsize=14, fontweight='bold')

    sm = plt.cm.ScalarMappable(cmap=plt.cm.RdBu, norm=plt.Normalize(vmin=-1,vmax=1))
    sm.set_array([])
    fig.colorbar(sm, ax=ax)

    plt.axis('off')
    plt.tight_layout()
    plt.savefig(save_path_svg, format='svg')
    if save_path_png:
        plt.savefig(save_path_png, dpi=200)
    plt.close()
    return save_path_png

# ==============================================
# GIF CREATION
# ==============================================
def create_gif(image_paths, gif_path):
    images = [imageio.imread(p) for p in image_paths]
    imageio.mimsave(gif_path, images, duration=1.2)

# ==============================================
# MAIN PIPELINE
# ==============================================
def main_pipeline(omega_df, traits_df, abundance_df, metadata_df, species_family_df):

    species_cols = omega_df['species_code'].tolist()
    print("Preprocessing data...")
    merged_df = preprocess_data(abundance_df, metadata_df, species_cols)

    builder = SectorGraphBuilder(omega_df, traits_df, merged_df, species_family_df)

    os.makedirs("outputs/svg", exist_ok=True)
    os.makedirs("outputs/png", exist_ok=True)
    os.makedirs("outputs/gif", exist_ok=True)
    os.makedirs("outputs/gephi", exist_ok=True)

    station_groups = merged_df.groupby('station')

    for station, df_station in station_groups:

        print(f"\nProcessing station: {station}")
        df_station = df_station.sort_values('year')
        image_paths = []

        for _, row in df_station.iterrows():
            name = row['station_year']
            svg_path = f"outputs/svg/{name}.svg"
            png_path = f"outputs/png/{name}.png"
            img_path = plot_network(builder, row, svg_path, png_path)
            image_paths.append(img_path)

            # EXPORT GEPHI
            G = nx.Graph()
            n = len(builder.species_names)
            G.add_nodes_from(builder.species_names)
            edge_index, edge_weights = builder.build_edge_index_and_weights()
            for i in range(edge_index.shape[1]):
                u = builder.species_names[edge_index[0,i]]
                v = builder.species_names[edge_index[1,i]]
                if u<v:
                    G.add_edge(u, v, weight=edge_weights[i].item())
            nx.write_gexf(G, f"outputs/gephi/{name}.gexf")

        # CREATE GIF
        gif_path = f"outputs/gif/{station}.gif"
        create_gif(image_paths, gif_path)
        print(f"✓ GIF created: {gif_path}")

# ==============================================
# RUN
# ==============================================
main_pipeline(omega_df, traits_df, abundance_df, metadata_df, species_family_df)
