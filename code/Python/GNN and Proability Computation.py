
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 18:05:37 2026

@author: p03565
"""

import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.cluster.hierarchy as sch

import pandas as pd
import gower
import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score
import numpy as np
#from sklearn_extra.cluster import KMedoids
import torch
import torch.nn.functional as F
from torch_geometric.data import Data
from torch_geometric.nn import GCNConv

# === 1. Load ===
# Replace with your file path
csv_file_path = 'C:\\Fish Data - Dr. Camara\\Moselle\\GTN_results_First\\NB LCMM\\GLMM - Baseline Stations -ZE1,2,3 and AFE 5,7\\Data Files\\Baseline_ZE1,2,3etAFE5,7_updated_11122025.csv'

df = pd.read_csv(csv_file_path)

df.columns

import pandas as pd
import numpy as np

# ---- identify species columns ----
# adjust this if you have a species name prefix
non_species_cols = [
'Date', 'Station', 'Sector', 'Campaign', 'Year', 'first_detection',
       'time_since_arrival', 'dale', 'flow_rate_low_water',
       'flow_rate_average_water', 'galet', 'gravier', 'helophytes_m_2_',
       'hydrophytes_m_2_', 'location_of_spaces__1to7__within_sectors',
       'median_height_of_the_seine_height', 'nearby_open_space', 'pente_',
       'pente_code', 'presence_nav_canal', 'presence_nav_channel', 'sable',
       'shore', 'space_qualification', 'substrate_m_2_',
       'surf_obs_hy_hl_ar_sNU_mÃƒ_Ã‚Â²_', 'total_height_of_the_seine',
       'trees_m_2_', 'vase', 'AverageDepth_m_', 'Bloc', 'Litere',
       'within_sector_sim', 'between_sector_sim',
       'distance_bw_station_downinfrastructure_km',
       'distance_bw_border_infra_downstream_sector',
       'distance_bw_border_infra_upstream_sector',
       'distance_bw_station_border', 'GTN'
]

species_cols = [c for c in df.columns if c not in non_species_cols]

species_cols
# ---- station-level richness ----
df['species_richness_no_GTN'] = (df[species_cols] > 0).sum(axis=1)

# inverse-distance accessibility
df['A_station_down'] = 1.0 / (1.0 + df['distance_bw_station_downinfrastructure_km'])

# scale distance to border (optional but recommended)
from sklearn.preprocessing import StandardScaler

scaler = StandardScaler()
df['P_station_border'] = scaler.fit_transform(
    df[['distance_bw_station_border']]
)

##Aggregate to sector-level
sector_agg = (
    df
    .groupby('Sector')
    .agg(
        mean_richness=('species_richness', 'mean'),
        mean_A_station_down=('A_station_down', 'mean'),
        mean_P_station_border=('P_station_border', 'mean')
    )
    .reset_index()
)

# scale distance to border
scaler = StandardScaler()
df['P_station_border'] = scaler.fit_transform(df[['distance_bw_station_border']])

# sector-level distances
df['A_border_down'] = 1.0 / (1.0 + df['distance_bw_border_infra_downstream_sector'])
df['A_border_up'] = 1.0 / (1.0 + df['distance_bw_border_infra_upstream_sector'])
df['A_directional_diff'] = df['A_border_down'] - df['A_border_up']
df['A_openness'] = df['A_border_down'] + df['A_border_up']

# GTN presence as binary
df['GTN_presence'] = (df['GTN'] < 0.5).astype(int)

# ----------------------------
# 4. Station-level GNN features (remove temporal/campaign info)
# ----------------------------
feature_cols = [
    'between_sector_sim',
    'within_sector_sim',
    'A_station_down',
    'P_station_border',
    'A_border_down',
    'A_border_up',
    'A_directional_diff',
    'A_openness',
    'GTN_presence'
]

X_all = df[feature_cols].copy()
scaler = StandardScaler()
X_all[feature_cols] = scaler.fit_transform(X_all[feature_cols])

X_tensor = torch.tensor(X_all.values, dtype=torch.float)
y_tensor = torch.tensor(df['species_richness_no_GTN'].values, dtype=torch.float)

# ----------------------------
# 5. Graph construction
# ----------------------------
# Node ID per station-year
df = df.reset_index(drop=True)
df['node_id'] = df.index
N = df.shape[0]

edge_index = []
edge_weight = []

# Temporal edges (optional, can be skipped since features are station-level)
for station, g in df.groupby('Station'):
    g = g.sort_values('Year')
    nodes = g['node_id'].values
    for i in range(len(nodes) - 1):
        n1, n2 = nodes[i], nodes[i + 1]
        edge_index.extend([[n1, n2], [n2, n1]])
        edge_weight.extend([1.0, 1.0])

# Spatial/similarity edges within sector
for (sector, year), g in df.groupby(['Sector', 'Year']):
    nodes = g['node_id'].values
    sims = g['within_sector_sim'].values
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            w = 0.5 * (sims[i] + sims[j])
            edge_index.extend([[nodes[i], nodes[j]], [nodes[j], nodes[i]]])
            edge_weight.extend([w, w])

edge_index = torch.tensor(edge_index, dtype=torch.long).t()
edge_weight = torch.tensor(edge_weight, dtype=torch.float)

data = Data(
    x=X_tensor,
    edge_index=edge_index,
    edge_attr=edge_weight,
    y=y_tensor
)

# ----------------------------
# 6. Define Weighted GCN
# ----------------------------
class StationLevelGCN(torch.nn.Module):
    def __init__(self, in_channels):
        super().__init__()
        self.conv1 = GCNConv(in_channels, 32)
        self.conv2 = GCNConv(32, 16)
        self.lin = torch.nn.Linear(16, 1)

    def forward(self, data):
        x, edge_index, edge_weight = data.x, data.edge_index, data.edge_attr
        x = F.relu(self.conv1(x, edge_index, edge_weight))
        x = F.relu(self.conv2(x, edge_index, edge_weight))
        out = self.lin(x)
        return out.squeeze(), x  # return embeddings

model = StationLevelGCN(in_channels=data.num_node_features)
optimizer = torch.optim.Adam(model.parameters(), lr=0.01, weight_decay=1e-4)
loss_fn = torch.nn.MSELoss()

# ----------------------------
# 7. Train GCN
# ----------------------------
n_epochs = 1000
for epoch in range(n_epochs):
    model.train()
    optimizer.zero_grad()
    preds, _ = model(data)
    loss = loss_fn(preds, data.y)
    loss.backward()
    optimizer.step()
    if epoch % 50 == 0:
        print(f"Epoch {epoch:03d} | Loss: {loss.item():.4f}")

# ----------------------------
# 8. Extract embeddings
# ----------------------------
model.eval()
with torch.no_grad():
    _, embeddings = model(data)
embeddings = embeddings.cpu().numpy()

# Attach to dataframe
emb_df = df[['Station', 'Year', 'Sector', 'GTN_presence', 'species_richness_no_GTN']].copy()
for i in range(embeddings.shape[1]):
    emb_df[f'emb_{i}'] = embeddings[:, i]
    
 # ----------------------------
# 9. UMAP for 2D visualization
# ----------------------------
import umap
reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=42)
emb_2d = reducer.fit_transform(embeddings)
emb_df['UMAP1'] = emb_2d[:,0]
emb_df['UMAP2'] = emb_2d[:,1]

# Plot by sector
plt.figure(figsize=(6,5))
for sector, group in emb_df.groupby('Sector'):
    plt.scatter(group['UMAP1'], group['UMAP2'], label=sector, s=30)
plt.xlabel('UMAP1')
plt.ylabel('UMAP2')
plt.title('GNN Embeddings by Sector (Station-level only)')
plt.legend()
plt.show()

# ----------------------------
# 10. Optional: Weighted degree
# ----------------------------
weighted_degree = np.zeros(len(df))
edges_to_plot = edge_index.cpu().numpy().T
edge_w = edge_weight.cpu().numpy()
for idx, (i,j) in enumerate(edges_to_plot):
    weighted_degree[i] += edge_w[idx]
    weighted_degree[j] += edge_w[idx]
emb_df['weighted_degree'] = weighted_degree

# ----------------------------
# 11. Network plot helper
# ----------------------------
def plot_network(df_plot, node_color, color_label, filename=None):
    plt.figure(figsize=(10,8))
    ax = plt.gca()
    for idx, (i,j) in enumerate(edges_to_plot):
        xi, yi = df_plot.loc[i, 'UMAP1'], df_plot.loc[i, 'UMAP2']
        xj, yj = df_plot.loc[j, 'UMAP1'], df_plot.loc[j, 'UMAP2']
        ax.plot([xi, xj], [yi, yj], color='gray', alpha=0.3, lw=0.5)
    scatter = ax.scatter(
        df_plot['UMAP1'], df_plot['UMAP2'],
        c=node_color,
        s=100*(df_plot['weighted_degree']/df_plot['weighted_degree'].max()+0.3),
        cmap='coolwarm',
        edgecolor='k', linewidth=0.3
    )
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label(color_label)
    plt.xlabel('UMAP1')
    plt.ylabel('UMAP2')
    plt.title(f'GNN Embeddings: {color_label}')
    plt.tight_layout()
    if filename:
        plt.savefig(filename, dpi=300)
    plt.show()

# Plot GTN presence
plot_network(emb_df, node_color=emb_df['GTN_presence'], color_label='GTN Presence', filename='GNN_GTNPresence_station.svg')

# Plot species richness
plot_network(emb_df, node_color=emb_df['species_richness_no_GTN'], color_label='Species Richness', filename='GNN_SpeciesRichness_station.svg')

#Embedded#
from sklearn.decomposition import PCA

def plot_network_clustered_with_variance(df_plot, embeddings, cluster_labels,
                                         node_color, color_label,
                                         top_n_labels=5, filename=None):
    """
    Plot GNN station-level network with selective node labels, cross-cluster edges,
    and PCA explained variance.

    Parameters
    ----------
    df_plot : pd.DataFrame
        Must contain 'UMAP1', 'UMAP2', 'weighted_degree', 'Station'.
    embeddings : np.array
        Original GNN embeddings (N_nodes x N_features) for PCA explained variance.
    cluster_labels : array-like
        Cluster assignment per node (e.g., ZE vs AFE).
    node_color : array-like
        Node color values.
    color_label : str
        Label for colorbar.
    top_n_labels : int, optional
        Number of nodes to label per cluster. Default=5.
    filename : str, optional
        Path to save figure.
    """
    # --- Compute PCA explained variance ---
    pca = PCA(n_components=2)
    pca_2d = pca.fit_transform(embeddings)
    explained_var_ratio = pca.explained_variance_ratio_
    total_explained = explained_var_ratio.sum()

    # --- Setup plot ---
    plt.figure(figsize=(10,8))
    ax = plt.gca()
    edges_to_plot_np = edge_index.cpu().numpy().T
    edge_w_np = edge_weight.cpu().numpy()

    # --- Draw edges ---
    for idx, (i, j) in enumerate(edges_to_plot_np):
        xi, yi = df_plot.loc[i, 'UMAP1'], df_plot.loc[i, 'UMAP2']
        xj, yj = df_plot.loc[j, 'UMAP1'], df_plot.loc[j, 'UMAP2']

        # Highlight cross-cluster edges
        if cluster_labels[i] != cluster_labels[j]:
            ax.plot([xi, xj], [yi, yj], color='purple', alpha=0.6, lw=1.2)
        else:
            ax.plot([xi, xj], [yi, yj], color='gray', alpha=0.2, lw=0.5)

    # --- Draw nodes ---
    node_sizes = 100 * (df_plot['weighted_degree'] / df_plot['weighted_degree'].max() + 0.3)
    scatter = ax.scatter(
        df_plot['UMAP1'], df_plot['UMAP2'],
        c=node_color,
        s=node_sizes,
        cmap='coolwarm',
        edgecolor='k', linewidth=0.3
    )

    # --- Colorbar ---
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label(color_label)

    # --- Label top nodes per cluster ---
    labeled_nodes = []
    for cluster in np.unique(cluster_labels):
        cluster_nodes = np.where(cluster_labels == cluster)[0]
        # rank by weighted degree
        top_nodes = cluster_nodes[np.argsort(df_plot.loc[cluster_nodes, 'weighted_degree'].values)[-top_n_labels:]]
        labeled_nodes.extend(top_nodes)

    for idx in labeled_nodes:
        xi, yi = df_plot.loc[idx, 'UMAP1'], df_plot.loc[idx, 'UMAP2']
        ax.text(xi, yi, df_plot.loc[idx, 'Station'],
                fontsize=8, fontweight='bold', alpha=0.9)

    plt.xlabel('UMAP1')
    plt.ylabel('UMAP2')
    plt.title(f'GNN Embeddings: {color_label}\nPCA 2D Explained Variance = {total_explained:.2f}')
    plt.tight_layout()

    if filename:
        plt.savefig(filename, dpi=300)
    plt.show()

from sklearn.cluster import KMeans

# 1. Cluster embeddings (2 main clusters: ZE vs AFE)
kmeans = KMeans(n_clusters=2, random_state=42).fit(embeddings)
emb_df['cluster'] = kmeans.labels_
cluster_labels = emb_df['cluster'].values

# 2. Plot GTN Presence with PCA variance info
plot_network_clustered_with_variance(
    df_plot=emb_df,
    embeddings=embeddings,
    cluster_labels=cluster_labels,
    node_color=emb_df['GTN_presence'],
    color_label='GTN Presence',
    top_n_labels=3,
    filename='GNN_GTNPresence_station_clustered_variance.svg'
)

# 3. Plot Species Richness similarly
plot_network_clustered_with_variance(
    df_plot=emb_df,
    embeddings=embeddings,
    cluster_labels=cluster_labels,
    node_color=emb_df['species_richness_no_GTN'],
    color_label='Species Richness',
    top_n_labels=3,
    filename='GNN_SpeciesRichness_station_clustered_variance.svg'
)

def plot_network_station_labels(df_plot, embeddings, cluster_labels,
                                node_color, color_label, filename=None):
    """
    Plot GNN station-level network with all station labels, cross-cluster edges,
    weighted node sizes, and PCA explained variance.

    Parameters
    ----------
    df_plot : pd.DataFrame
        Must contain 'UMAP1', 'UMAP2', 'weighted_degree', 'Station'.
    embeddings : np.array
        Original GNN embeddings (N_nodes x N_features) for PCA explained variance.
    cluster_labels : array-like
        Cluster assignment per node (e.g., ZE vs AFE).
    node_color : array-like
        Node color values.
    color_label : str
        Label for colorbar.
    filename : str, optional
        Path to save figure.
    """
    # --- Compute PCA explained variance ---
    pca = PCA(n_components=2)
    pca.fit(embeddings)
    total_explained = pca.explained_variance_ratio_.sum()

    # --- Setup plot ---
    plt.figure(figsize=(12,10))
    ax = plt.gca()
    edges_to_plot_np = edge_index.cpu().numpy().T
    edge_w_np = edge_weight.cpu().numpy()

    # --- Draw edges ---
    for idx, (i, j) in enumerate(edges_to_plot_np):
        xi, yi = df_plot.loc[i, 'UMAP1'], df_plot.loc[i, 'UMAP2']
        xj, yj = df_plot.loc[j, 'UMAP1'], df_plot.loc[j, 'UMAP2']

        # Highlight cross-cluster edges
        if cluster_labels[i] != cluster_labels[j]:
            ax.plot([xi, xj], [yi, yj], color='purple', alpha=0.6, lw=1.2)
        else:
            ax.plot([xi, xj], [yi, yj], color='gray', alpha=0.2, lw=0.5)

    # --- Draw nodes ---
    node_sizes = 100 * (df_plot['weighted_degree'] / df_plot['weighted_degree'].max() + 0.3)
    scatter = ax.scatter(
        df_plot['UMAP1'], df_plot['UMAP2'],
        c=node_color,
        s=node_sizes,
        cmap='coolwarm',
        edgecolor='k', linewidth=0.3
    )

    # --- Colorbar ---
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label(color_label)

    # --- Label every node with its Station ---
    for idx in range(df_plot.shape[0]):
        xi, yi = df_plot.loc[idx, 'UMAP1'], df_plot.loc[idx, 'UMAP2']
        ax.text(xi, yi, df_plot.loc[idx, 'Station'],
                fontsize=8, fontweight='bold', alpha=0.9)

    plt.xlabel('UMAP1')
    plt.ylabel('UMAP2')
    plt.title(f'GNN Embeddings: {color_label}\nPCA 2D Explained Variance = {total_explained:.2f}')
    plt.tight_layout()

    if filename:
        plt.savefig(filename, dpi=300)
    plt.show()

# Assuming embeddings and clusters are already computed
plot_network_station_labels(
    df_plot=emb_df,
    embeddings=embeddings,
    cluster_labels=emb_df['cluster'].values,
    node_color=emb_df['GTN_presence'],
    color_label='GTN Presence',
    filename='GNN_GTNPresence_station_labeled.svg'
)

plot_network_station_labels(
    df_plot=emb_df,
    embeddings=embeddings,
    cluster_labels=emb_df['cluster'].values,
    node_color=emb_df['species_richness_no_GTN'],
    color_label='Species Richness',
    filename='GNN_SpeciesRichness_station_labeled.svg'
)

# -------------------------------
# 1️⃣ Recompute GTN Spread Probability
# -------------------------------

# Distance-based connectivity
# Use previously computed variables: A_station_down, P_station_border, A_border_down, A_border_up
# Normalize to [0,1]
distance_factors = ['A_station_down', 'P_station_border', 'A_border_down', 'A_border_up']
for col in distance_factors:
    df[col+'_norm'] = (df[col] - df[col].min()) / (df[col].max() - df[col].min())

# Trait-based spread factor
# Assign weights for traits (tunable)
trait_weights = {
    'Benthic': 0.7,       # Benthic habitat reduces movement compared to pelagic
    'LMAX': 1.0,          # smaller fish = higher probability
    'Migration': 0.5,     # Resident = lower dispersal across segments
    'Rheophily': 0.9      # EUR = generally tolerant, higher probability
}

# Compute trait factor for GTN (all stations)
# normalized between 0 and 1
trait_factor = (
    trait_weights['Benthic'] *
    trait_weights['LMAX'] *
    trait_weights['Migration'] *
    trait_weights['Rheophily']
)
# scale to [0,1]
trait_factor = trait_factor / max(trait_weights.values())**len(trait_weights)

# Combine distance and trait factors for each station → probability of spread
df['GTN_spread_prob'] = (
    df[['A_station_down_norm', 'P_station_border_norm', 'A_border_down_norm', 'A_border_up_norm']].mean(axis=1)
) * trait_factor

# Ensure probability is between 0 and 1
df['GTN_spread_prob'] = df['GTN_spread_prob'].clip(0,1)

# Quick sanity check
print(df[['Station','GTN_spread_prob']])

# -------------------------------
# 2️⃣ Prepare GNN embeddings and save
# -------------------------------
# -------------------------------
# -------------------------------
model.eval()
with torch.no_grad():
    out, embeddings = model(data)

emb_df = df[['Station', 'Campaign', 'Year']].copy()
emb_cols = [f'emb_{i}' for i in range(embeddings.shape[1])]
emb_df[emb_cols] = embeddings.cpu().numpy()

emb_df.to_csv('GNN_embeddings_station_level.csv', index=False)
print("GNN embeddings saved to GNN_embeddings_station_level.csv")
# -------------------------------
# 3️⃣ Optional: Visualize GTN Spread Probabilities on UMAP
# -------------------------------
emb_df['GTN_spread_prob'] = df['GTN_spread_prob'].values

import umap
reducer = umap.UMAP(n_neighbors=5, min_dist=0.3, random_state=42)
emb_2d = reducer.fit_transform(embeddings.cpu().numpy())

emb_df['UMAP1'] = emb_2d[:,0]
emb_df['UMAP2'] = emb_2d[:,1]

plt.figure(figsize=(12,8))
scatter = plt.scatter(
    emb_df['UMAP1'], emb_df['UMAP2'],
    c=emb_df['GTN_spread_prob'],
    cmap='Reds', s=150, edgecolor='k', alpha=0.9
)
plt.colorbar(scatter, label='GTN Spread Probability')
for i, row in emb_df.iterrows():
    plt.text(row['UMAP1'], row['UMAP2'], row['Station'], fontsize=9, ha='center', va='bottom')
plt.xlabel('UMAP1')
plt.ylabel('UMAP2')
plt.title('Station-level GNN Embeddings with Recomputed GTN Spread Probability')
plt.tight_layout()
plt.show()

import matplotlib.pyplot as plt
import numpy as np

plt.figure(figsize=(12, 8))

# Scatter all nodes
plt.scatter(
    emb_df['UMAP1'],
    emb_df['UMAP2'],
    c=emb_df['GTN_spread_prob'],
    cmap='Reds',
    s=60,
    alpha=0.6,
    edgecolor='k'
)

# ---- Label station centroids ----
for station, g in emb_df.groupby('Station'):
    cx = g['UMAP1'].mean()
    cy = g['UMAP2'].mean()
    plt.text(
        cx, cy, station,
        fontsize=12,
        fontweight='bold',
        ha='center', va='center',
        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="black", alpha=0.8)
    )

plt.colorbar(label='GTN Spread Probability')
plt.xlabel('UMAP1')
plt.ylabel('UMAP2')
plt.title('UMAP of GNN Station Embeddings (Station-level labels)')
plt.tight_layout()
plt.show()


#Add edges
river_order_map = {
    "ZE1": 1,
    "ZE2": 2,
    "ZE3": 3,
    "AFE5": 4,
    "AFE7": 5
}

emb_df["river_order"] = emb_df["Station"].map(river_order_map)

station_centroids = (
    emb_df
    .groupby("Station")[["UMAP1", "UMAP2", "river_order"]]
    .mean()
    .reset_index()
    .sort_values("river_order")
)

import matplotlib.pyplot as plt

plt.figure(figsize=(9, 7))

# ---- all points (faint) ----
plt.scatter(
    emb_df["UMAP1"],
    emb_df["UMAP2"],
    s=35,
    alpha=0.4
)

# ---- station centroids ----
plt.scatter(
    station_centroids["UMAP1"],
    station_centroids["UMAP2"],
    s=180,
    zorder=3
)

# ---- label every station (now readable) ----
for _, row in station_centroids.iterrows():
    plt.text(
        row["UMAP1"],
        row["UMAP2"],
        row["Station"],
        fontsize=10,
        fontweight="bold",
        ha="center",
        va="center"
    )

# ---- arrows: downstream → upstream ----
for i in range(len(station_centroids) - 1):
    x1, y1 = station_centroids.iloc[i][["UMAP1", "UMAP2"]]
    x2, y2 = station_centroids.iloc[i + 1][["UMAP1", "UMAP2"]]

    plt.arrow(
        x1, y1,
        x2 - x1,
        y2 - y1,
        length_includes_head=True,
        head_width=0.03,
        head_length=0.05,
        linewidth=2
    )

plt.xlabel("UMAP-1")
plt.ylabel("UMAP-2")
plt.title("Directional GTN spread along river continuum")
plt.tight_layout()
plt.show()

#Another
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

emb_df.columns 
station_centroids = (
    emb_df
    .groupby("Station")[["UMAP1", "UMAP2", "GTN_spread_prob"]]
    .mean()
    .reset_index()
)

river_order = ["ZE1", "ZE2", "ZE3", "AFE5", "AFE7"]
station_centroids["order"] = station_centroids["Station"].apply(
    lambda x: river_order.index(x)
)
station_centroids = station_centroids.sort_values("order")

plt.figure(figsize=(9, 7))

# ---- background points (keep same UMAP layout) ----
plt.scatter(
    emb_df["UMAP1"],
    emb_df["UMAP2"],
    s=30,
    alpha=0.35
)

# ---- station centroids ----
plt.scatter(
    station_centroids["UMAP1"],
    station_centroids["UMAP2"],
    s=220,
    zorder=3
)

# ---- station labels + GTN probability ----
for _, row in station_centroids.iterrows():
    label = f"{row['Station']}\nP={row['GTN_spread_prob']:.2f}"
    plt.text(
        row["UMAP1"],
        row["UMAP2"],
        label,
        fontsize=10,
        fontweight="bold",
        ha="center",
        va="center"
    )

# ---- directional arrows (downstream → upstream) ----
for i in range(len(station_centroids) - 1):
    start = station_centroids.iloc[i]
    end = station_centroids.iloc[i + 1]

    arrow = FancyArrowPatch(
        (start["UMAP1"], start["UMAP2"]),
        (end["UMAP1"], end["UMAP2"]),
        arrowstyle="->",
        linewidth=2,
        mutation_scale=18,
        zorder=4
    )
    plt.gca().add_patch(arrow)

plt.xlabel("UMAP-1")
plt.ylabel("UMAP-2")
plt.title("GTN spread probability and directional propagation along river network")
plt.tight_layout()
plt.show()

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
import matplotlib as mpl

# ---- colormap for GTN probability ----
cmap = plt.cm.Reds
norm = mpl.colors.Normalize(
    vmin=station_centroids["GTN_spread_prob"].min(),
    vmax=station_centroids["GTN_spread_prob"].max()
)

plt.figure(figsize=(9, 7))

# ---- background points (same layout, faint) ----
plt.scatter(
    emb_df["UMAP1"],
    emb_df["UMAP2"],
    s=25,
    alpha=0.25,
    color="grey"
)

# ---- station centroids (colored by probability) ----
sc = plt.scatter(
    station_centroids["UMAP1"],
    station_centroids["UMAP2"],
    c=station_centroids["GTN_spread_prob"],
    cmap=cmap,
    norm=norm,
    s=260,
    edgecolor="black",
    zorder=3
)

# ---- labels: station + probability ----
for _, row in station_centroids.iterrows():
    label = f"{row['Station']}\nP={row['GTN_spread_prob']:.2f}"
    plt.text(
        row["UMAP1"],
        row["UMAP2"],
        label,
        fontsize=10,
        fontweight="bold",
        ha="center",
        va="center",
        zorder=4
    )

# ---- directional arrows (downstream → upstream) ----
for i in range(len(station_centroids) - 1):
    start = station_centroids.iloc[i]
    end = station_centroids.iloc[i + 1]

    arrow = FancyArrowPatch(
        (start["UMAP1"], start["UMAP2"]),
        (end["UMAP1"], end["UMAP2"]),
        arrowstyle="->",
        linewidth=2.2,
        mutation_scale=18,
        color="black",
        zorder=4
    )
    plt.gca().add_patch(arrow)

# ---- colorbar ----
cbar = plt.colorbar(sc, pad=0.02)
cbar.set_label("GTN spread probability", fontsize=11)

plt.xlabel("UMAP-1")
plt.ylabel("UMAP-2")
plt.title("Directional GTN spread probability along the river network")
plt.tight_layout()
plt.show()

#Final plot hopefully
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

plt.figure(figsize=(14, 10))

# ---- Scatter all nodes (same layout, slightly cleaner) ----
plt.scatter(
    emb_df['UMAP1'],
    emb_df['UMAP2'],
    c=emb_df['GTN_spread_prob'],
    cmap='Reds',
    s=45,               # slightly smaller
    alpha=0.5,
    edgecolor='k',
    linewidth=0.3
)

# ---- Compute station centroids ----
centroids = (
    emb_df
    .groupby('Station')[['UMAP1', 'UMAP2', 'GTN_spread_prob']]
    .mean()
    .reset_index()
)

# ---- River order (explicit, not inferred) ----
river_order = ['ZE1', 'ZE2', 'ZE3', 'AFE5', 'AFE7']
centroids['order'] = centroids['Station'].apply(
    lambda x: river_order.index(x)
)
centroids = centroids.sort_values('order')

# ---- Label station centroids WITH probability ----
for _, row in centroids.iterrows():
    label = f"{row['Station']}\nP={row['GTN_spread_prob']:.2f}"
    plt.text(
        row['UMAP1'],
        row['UMAP2'],
        label,
        fontsize=11,
        fontweight='bold',
        ha='center',
        va='center',
        bbox=dict(
            boxstyle="round,pad=0.25",
            fc="white",
            ec="black",
            alpha=0.85
        ),
        zorder=4
    )

# ---- Directional arrows (downstream → upstream) ----
for i in range(len(centroids) - 1):
    start = centroids.iloc[i]
    end = centroids.iloc[i + 1]

    arrow = FancyArrowPatch(
        (start['UMAP1'], start['UMAP2']),
        (end['UMAP1'], end['UMAP2']),
        arrowstyle='->',
        mutation_scale=20,   # makes it pointy
        linewidth=2.2,
        color='black',
        zorder=3
    )
    plt.gca().add_patch(arrow)

# ---- Colorbar and labels ----
cbar = plt.colorbar(label='GTN Spread Probability')
cbar.ax.tick_params(labelsize=10)

plt.xlabel('UMAP1')
plt.ylabel('UMAP2')
plt.title('UMAP of GNN Station Embeddings with GTN Spread Probability')
plt.tight_layout()
plt.show()

#FI Plots#
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from copy import deepcopy

# ---- Step 1: baseline predictions ----
model.eval()
with torch.no_grad():
    base_preds, _ = model(data)  # assumes StationLevelGCN forward returns (out, embeddings)
    base_preds = base_preds.cpu().numpy()

# ---- Step 2: define features to test ----
# Only distance & GTN traits, exclude Year/Campaign
feat_cols = [
    'between_sector_sim', 'within_sector_sim',
    'A_station_down', 'P_station_border',
    'A_border_down', 'A_border_up',
    'A_directional_diff', 'A_openness',
]

# If your X_tensor columns are in the same order as X_all.columns
feat_indices = [X_all.columns.get_loc(f) for f in feat_cols]

# ---- Step 3: permutation importance ----
drop_importances = {}

for idx, col in zip(feat_indices, feat_cols):
    X_perm = X_tensor.clone()
    # Shuffle the column
    shuffled = X_perm[:, idx][torch.randperm(X_perm.size(0))]
    X_perm[:, idx] = shuffled
    # Recompute predictions
    with torch.no_grad():
        preds_perm, _ = model(Data(x=X_perm, edge_index=data.edge_index, edge_attr=data.edge_attr))
        preds_perm = preds_perm.cpu().numpy()
    # Mean absolute change
    delta = np.mean(np.abs(preds_perm - base_preds))
    drop_importances[col] = delta
    print(f"{col}: Δ={delta:.4f}")

# ---- Step 4: convert to DataFrame for plotting ----
importances_df = pd.DataFrame.from_dict(drop_importances, orient='index', columns=['MAE_drop'])
importances_df = importances_df.sort_values('MAE_drop', ascending=True)

# ---- Step 5: plot ----
plt.figure(figsize=(10,6))
plt.barh(importances_df.index, importances_df['MAE_drop'], color='teal')
plt.xlabel('Mean Absolute Change in GTN Spread Probability')
plt.title('Permutation Feature Importance (Distance & GTN Traits)')
plt.tight_layout()
plt.show()

#Save the emB_df
emb_df.to_csv("GNN_embeddings_with_Probaility_GTN_SPREAD.csv", index = False)
