import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import scipy.sparse
from matplotlib.colors import hex2color, ListedColormap
from matplotlib.patches import Patch
import anndata
from scipy.cluster.hierarchy import leaves_list
from sklearn.preprocessing import OrdinalEncoder
from scipy.cluster.hierarchy import dendrogram
import sys
import os
from numpy_groupies.aggregate_numpy import aggregate
import cytograph as cg
import cytograph.visualization as cgplot
#import cytograph.plotting as cgplot
import loompy
adata=anndata.read_h5ad("pig.h5ad")
#################################################
# Helper functions
def sparkline(ax, x, ymax, color, plot_label):
    n_clusters = x.shape[0]
    if ymax is None:
        ymax = np.max(x)
    ax.bar(np.arange(n_clusters) + 0.5, x, color=color, width=1, lw=0)
    ax.set_xlim(0, n_clusters + 1)
    ax.set_ylim(0, ymax)
    ax.axis("off")
    ax.text(0, 0, plot_label, va="bottom", ha="right", transform=ax.transAxes)



def plot_genes(ax, mgenes, mean_x, genes):
    # Add the markers
    m = []
    m_names = []
    for gene in mgenes:
        gene_ix = np.where(genes == gene)[0][0]
        m.append(mean_x[:, gene_ix])
        m_names.append(f"{gene}")
    n_genes = len(m_names)
    n_clusters = labels.max() + 1
    # Normalize
    x = np.array(m)
    x_norm = cg.div0(x.T, np.percentile(x, 99.9, axis=1)).T
    bg = np.zeros_like(x_norm) + 0.9
    x_norm = np.ma.masked_where(x == 0, x_norm)
    ax.imshow(bg, vmin=0, vmax=1, cmap=plt.cm.gray, aspect="auto", extent=(0, n_clusters, n_genes, 0), alpha=1, interpolation="nearest", resample=False)
    ax.imshow(x_norm, cmap="inferno_r", vmax=1, interpolation="nearest", aspect="auto", alpha=1, extent=(0, n_clusters, n_genes, 0), resample=False)
    #ax.imshow(np.log10(x_norm + 0.001), vmin=-1, vmax=2, cmap="RdGy_r", interpolation="none", aspect="auto", extent=(0, n_clusters, n_genes, 0))
    ax.set_yticks(np.arange(len(m_names)) + 0.5)
    ax.set_yticklabels(m_names, fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    
import numpy as np
import pandas as pd
from sklearn.preprocessing import OrdinalEncoder


adata.obs['subclass_encoded'] = adata.obs['subclass'].astype('category').cat.codes
unique_subclasses = adata.obs['subclass'].astype('category').cat.categories
n_clusters = len(unique_subclasses)


rois = adata.obs['region'].dropna().unique() 
n_rois = len(rois)
le = OrdinalEncoder(categories=[rois])
le = OrdinalEncoder(categories=[['PFC', 'STR', 'THA']])


le.fit(adata.obs[['region']])  


roidistro = np.zeros((n_rois, n_clusters), dtype=int)


for label in range(n_clusters):

    mask = (adata.obs['subclass_encoded'] == label)

    if mask.sum() == 0:
        print(f"Warning: Cluster {unique_subclasses[label]} has no cells.")
        continue

    regions = adata.obs.loc[mask, 'region'].dropna().values
    if len(regions) == 0:
        print(f"Warning: Cluster {unique_subclasses[label]} has no region info.")
        continue
    encoded = le.transform(regions.reshape(-1, 1)).flatten().astype(int) 
    counts = np.bincount(encoded, minlength=n_rois)
    roidistro[:, label] = counts



####################
import numpy as np
import pandas as pd
donors = ['L-GE-14-2', 'L-GE-14-4', 'L-GE-15-1', 'L-GE-15-2', 'L-LGE-13-1', 'L-LGE-13-2', 'L-STR-11-1', 'L-Str-4-1', 'L-Str-4-3', 'L-Str-5-2', 'L-Str-5-7', 'L-Str-6-2', 'L-Str-7-1', 'L-Str-7-2', 'L-Str-8-2', 'L-THA-11-1', 'L-THA-11-3', 'L-THA-13-2', 'L-THA-14-2', 'L-THA-14-4', 'L-THA-15-1', 'L-THA-15-2', 'L-Tha-4-1', 'L-Tha-4-3', 'L-Tha-5-2', 'L-Tha-5-7', 'L-Tha-6-1', 'L-Tha-6-2', 'L-Tha-7-1', 'L-Tha-7-2', 'L-Tha-8-1', 'L-Tha-8-2', 'L_PFC_4_1', 'L_PFC_4_3', 'L_PFC_5_2', 'L_PFC_5_7', 'L_PFC_6_1', 'L_PFC_6_2', 'L_PFC_7_1', 'L_PFC_7_2', 'PFC_9_1', 'PFC_11_1', 'PFC_11_3', 'PFC_13_1', 'PFC_13_2', 'PFC_14_2', 'PFC_14_4', 'PFC_15_1', 'PFC_15_2', 'R-STR-11-3']
donordist = np.zeros((n_clusters, len(donors)), dtype="float32")

for ix, donor in enumerate(donors):
    donor_cells = adata.obs['orig_ident'] == donor
    
    subclass_codes = pd.Categorical(adata.obs['subclass'][donor_cells]).codes
    
    donor_counts = np.bincount(subclass_codes, minlength=n_clusters)
    
    donordist[:, ix] = donor_counts


donordist = (donordist.T / donordist.sum(axis=1)).T
###########

times = adata.obs["time"].dropna().unique()
n_times = len(times)
le_time = OrdinalEncoder(categories=[times])
le_time.fit(adata.obs[["time"]])

time_dist = np.zeros((n_clusters, n_times), dtype="float32")
for label in range(n_clusters):
    mask = (adata.obs["subclass_encoded"] == label)
    time_values = adata.obs.loc[mask, "time"].dropna().values
    if len(time_values) == 0:
        continue
    encoded = le_time.transform(time_values.reshape(-1, 1)).flatten().astype(int)
    counts = np.bincount(encoded, minlength=n_times)
    time_dist[label, :] = counts

time_order = ['E45', 'E55', 'E66', 'E76', 'E85', 'E94', 'E104', 'E109', 'P0', 'P3']
n_times = len(time_order)

le_time = OrdinalEncoder(categories=[time_order])
le_time.fit(adata.obs[["time"]])

time_dist = np.zeros((n_clusters, n_times), dtype="float32")
for label in range(n_clusters):
    mask = (adata.obs["subclass_encoded"] == label)
    time_values = adata.obs.loc[mask, "time"].dropna().values
    if len(time_values) == 0:
        continue
    encoded = le_time.transform(time_values.reshape(-1, 1)).flatten().astype(int)
    counts = np.bincount(encoded, minlength=n_times)
    time_dist[label, :] = counts
    

time_dist = (time_dist.T / time_dist.sum(axis=1)).T
time_dist = np.nan_to_num(time_dist, nan=0)

################
genes = adata.var_names.values  

import scanpy as sc
sc.pp.log1p(adata) 
sc.tl.rank_genes_groups(adata, groupby='subclass', use_raw=False)  
mean_x = np.zeros((len(adata.obs['subclass'].cat.categories), adata.n_vars))
for i, cl in enumerate(adata.obs['subclass'].cat.categories):
    mean_x[i] = adata[adata.obs['subclass'] == cl].X.mean(axis=0)

n_cells_per_cluster = adata.obs['subclass'].value_counts().sort_index().values


labels = adata.obs['subclass'].cat.codes.values


rg = adata.obs['region'].copy()

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['svg.fonttype'] = 'none'

subclass_color_dict = {
    "RG": "#D62728",
    "RG_neuron": "#DB3A3B",
    "enIPC": "#E04E4E",
    "enIPC cycle": "#E56162",
    "inIPC": "#EA7575",
    "inIPC cycle": "#EF8888",
    "gIPC": "#F49B9C",
    "gIPC cycle": "#F9AFAF",
    "NPCs_CHGB": "#FFC3C3",
    "mExN2": "#1F77B4",
    "mExN1": "#277CB8",
    "DpL ExN": "#3082BC",
    "UpL ExN": "#3987C1",
    "L2_3 IT": "#428DC5",
    "L3_5 IT": "#4B92CA",
    "L5 ET": "#5398CE",
    "L5_6 IT": "#5C9DD2",
    "L5_6 NP": "#65A3D7",
    "L6 B": "#6EA8DB",
    "L6 CT": "#77AEE0",
    "ExN_Str": "#80B3E4",
    "ExN_CTX_FOXP2": "#88B9E8",
    "ExN_GPR39": "#91BEED",
    "ExN_GRIK1": "#9AC4F1",
    "ExN_HB": "#A3C9F6",
    "ExN_TAFA4": "#ACCFFA",
    "ExN_VS": "#B5D5FF",
    "cgeInN": "#2CA02C",
    "lgeInN": "#3CAD3C",
    "mgeInN": "#4CBB4C",
    "InN_GRIK1": "#5CC95C",
    "InN_MB": "#6CD66C",
    "InN_P3": "#7CE47C",
    "InN_TEL": "#8CF28C",     
    "dSPN": "#9467BD",
    "iSPN": "#BC99DE",
    "eSPN": "#E4CCFF",
    "Astro": "#FF7F0E",
    "OPC": "#FF8E2E",
    "OPC cycle": "#FF9D4E",
    "COP": "#FFAD6F",
    "Oligo": "#FFBC8F",
    "Micro": "#FFCBAF",
    "Micro cycle": "#FFDBD0",
    "CR": "#8C564B",
    "Ependymal": "#A77066",
    "END": "#C38B81",
    "Mural": "#DEA59C",
    "VLMC": "#FAC0B7",
    "unknown": "#BBBBBB"
}

class_order = ["PC", "ExN", "InN", "SPN", "Glia", "Others", "unknown"]
class_color_map = {
    "PC": "#d62728",
    "ExN": "#1f77b4",
    "InN": "#2ca02c",
    "SPN": "#9467bd",
    "Glia": "#ff7f0e",
    "Others": "#8c564b",
    "unknown": "#9E9E9E"
}

existing_subclasses = [s for s in subclass_color_dict.keys() if s in adata.obs['subclass'].cat.categories]
order_indices = [list(adata.obs['subclass'].cat.categories).index(s) for s in existing_subclasses]

if 'class' not in adata.obs:
    raise KeyError("adata.obs does not contain 'class' column. Please add class information.")
    
subclass_to_class = {}
for subclass in existing_subclasses:
    class_value = adata.obs[adata.obs['subclass'] == subclass]['class'].iloc[0]
    subclass_to_class[subclass] = class_value


class_order_list = [subclass_to_class[cl] for cl in existing_subclasses]

if 'total_counts' not in adata.obs:
    if scipy.sparse.issparse(adata.X):
        adata.obs['total_counts'] = adata.X.sum(axis=1).A1
    else:
        adata.obs['total_counts'] = adata.X.sum(axis=1)

if 'n_genes' not in adata.obs:
    if scipy.sparse.issparse(adata.X):
        adata.obs['n_genes'] = (adata.X > 0).sum(axis=1).A1
    else:
        adata.obs['n_genes'] = (adata.X > 0).sum(axis=1)

n_cells_per_cluster = adata.obs['subclass'].value_counts().sort_index().values
avg = np.zeros(len(existing_subclasses))
avg_genes = np.zeros(len(existing_subclasses))

for i, cl in enumerate(existing_subclasses):
    mask = adata.obs['subclass'] == cl
    avg[i] = np.mean(adata.obs.loc[mask, 'total_counts'])
    avg_genes[i] = np.mean(adata.obs.loc[mask, 'n_genes'])

# 4. 计算性别比例
# 确保存在sex列
if 'sex' in adata.obs:
    sex_categories = adata.obs['sex'].dropna().unique()
    n_sex = len(sex_categories)
    sex_dist = np.zeros((len(existing_subclasses), n_sex), dtype="float32")
    
    le_sex = OrdinalEncoder(categories=[sex_categories])
    le_sex.fit(adata.obs[["sex"]])
    
    for j, cl in enumerate(existing_subclasses):
        mask = (adata.obs['subclass'] == cl)
        sex_values = adata.obs.loc[mask, "sex"].dropna().values
        if len(sex_values) > 0:
            encoded = le_sex.transform(sex_values.reshape(-1, 1)).flatten().astype(int)
            counts = np.bincount(encoded, minlength=n_sex)
            sex_dist[j, :] = counts
    sex_dist = (sex_dist.T / sex_dist.sum(axis=1)).T
    sex_dist = np.nan_to_num(sex_dist, nan=0)
    sex_colors = []
    if 'Male' in sex_categories and 'Female' in sex_categories:
        sex_colors = ["#1f77b4", "#ff6b6b"]  
    else:
        sex_colors = plt.cm.tab10.colors[:n_sex]
else:
    print("Warning: 'sex' column not found in adata.obs. Skipping sex ratio plot.")
    sex_dist = None


n_cells_ordered = n_cells_per_cluster[order_indices]


rois = ['PFC', 'STR', 'THA']
le = OrdinalEncoder(categories=[rois])
le.fit(adata.obs[['region']])
roidistro = np.zeros((len(rois), len(existing_subclasses)), dtype=int)

for j, cl in enumerate(existing_subclasses):
    mask = (adata.obs['subclass'] == cl)
    regions = adata.obs.loc[mask, 'region'].dropna().values
    if len(regions) > 0:
        encoded = le.transform(regions.reshape(-1, 1)).flatten().astype(int)
        counts = np.bincount(encoded, minlength=len(rois))
        roidistro[:, j] = counts


time_order_original = ['E45', 'E55', 'E66', 'E76', 'E85', 'E94', 'E104', 'E109', 'P0', 'P3']
time_order = ['E45', 'E55', 'E66', 'E76', 'E85', 'E94', 'E104', 'E109', 'P0/P3']


le_time = OrdinalEncoder(categories=[time_order_original])
le_time.fit(adata.obs[["time"]])
time_dist = np.zeros((len(existing_subclasses), len(time_order_original)), dtype="float32")

for j, cl in enumerate(existing_subclasses):
    mask = (adata.obs['subclass'] == cl)
    time_values = adata.obs.loc[mask, "time"].dropna().values
    if len(time_values) > 0:
        encoded = le_time.transform(time_values.reshape(-1, 1)).flatten().astype(int)
        counts = np.bincount(encoded, minlength=len(time_order_original))
        time_dist[j, :] = counts


p0_p3 = time_dist[:, -2:].sum(axis=1, keepdims=True)
time_dist = np.hstack([time_dist[:, :-2], p0_p3])


time_dist = (time_dist.T / time_dist.sum(axis=1)).T
time_dist = np.nan_to_num(time_dist, nan=0)
time_dist_ordered = time_dist


sc.pp.log1p(adata)
mean_x = np.zeros((len(existing_subclasses), adata.n_vars))
for i, cl in enumerate(existing_subclasses):
    mean_x[i] = adata[adata.obs['subclass'] == cl].X.mean(axis=0)
genes = adata.var_names.values

n_rows = 8
height_ratios = [0.1, 0.1, 0.1, 0.1, 0.15, 0.08, 2.5, 0.5] 

fig, axes = plt.subplots(
    nrows=n_rows, 
    ncols=1, 
    sharex=True, 
    gridspec_kw={
        "height_ratios": height_ratios,
        "hspace": 0.05
    },
    figsize=(15, 8 + (1 if sex_dist is not None else 0)) 
)

def sparkline(ax, x, ymax, color, plot_label):
    n = x.shape[0]
    if ymax is None:
        ymax = np.max(x) * 1.1  
    ax.bar(np.arange(n) + 0.5, x, color=color, width=1, lw=0)
    ax.set_xlim(0, n + 1)
    ax.set_ylim(0, ymax)
    ax.axis("off")
    ax.spines['bottom'].set_visible(False)
    ax.text(0, 0, plot_label, va="bottom", ha="right", transform=ax.transAxes)

ax = axes[0]
class_colors = [class_color_map[cls] for cls in class_order_list]

for i, color in enumerate(class_colors):
    ax.add_patch(plt.Rectangle((i, 0), 1, 1, color=color))
    
    

ax.set_xlim(0, len(existing_subclasses))
ax.set_ylim(0, 1)
ax.axis("off")

ax.spines['bottom'].set_visible(False)

legend_handles = []
for cls, color in class_color_map.items():
    legend_handles.append(Patch(color=color, label=cls))
    

ax.legend(handles=legend_handles, loc='upper center', bbox_to_anchor=(0.5, 1.5), 
          ncol=len(class_color_map), fontsize=8, frameon=False)

ax = axes[1]
subclass_colors = [subclass_color_dict[cl] for cl in existing_subclasses]

for i, color in enumerate(subclass_colors):
    ax.add_patch(plt.Rectangle((i, 0), 1, 1, color=color))
    



ax.set_xlim(0, len(existing_subclasses))
ax.set_ylim(0, 1)
ax.axis("off")
ax.spines['bottom'].set_visible(False)
ax.text(0.01, 0.5, "Subclass Colors", va="center", ha="left", transform=ax.transAxes, fontsize=8)
sparkline(axes[2], n_cells_ordered, None, "orange", "# cells")
sparkline(axes[3], avg_genes, None, "green", "avg genes")
ax = axes[4]
roi_colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]
roidistro_norm = roidistro / np.maximum(roidistro.sum(axis=0), 1)  

color_matrix = np.zeros((len(rois), len(existing_subclasses), 4))
for i in range(len(rois)):
    rgb_color = hex2color(roi_colors[i])
    for j in range(len(existing_subclasses)):
        color_matrix[i, j, :3] = rgb_color
        color_matrix[i, j, 3] = roidistro_norm[i, j]




ax.imshow(color_matrix, aspect='auto', interpolation='nearest', origin='upper',
          extent=(0, len(existing_subclasses), len(rois), 0))
ax.set_yticks(np.arange(len(rois)) + 0.5)
ax.set_yticklabels(rois, fontsize=8)

for spine in ['top', 'right', 'bottom', 'left']:
    ax.spines[spine].set_visible(False)


if sex_dist is not None:
    ax = axes[5]
    bottom = np.zeros(len(existing_subclasses))
    for i in range(n_sex):
        ax.bar(
            np.arange(len(existing_subclasses)) + 0.5,
            sex_dist[:, i],
            width=1,
            bottom=bottom,
            color=sex_colors[i],
            edgecolor='none'
        )
        bottom += sex_dist[:, i]
    
    ax.set_xlim(0, len(existing_subclasses))
    ax.set_ylim(0, 1)
    ax.set_ylabel('Sex ratio')
    for spine in ['top', 'right', 'bottom', 'left']:
        ax.spines[spine].set_visible(False)
    
    legend_handles = []
    for i, sex in enumerate(sex_categories):
        legend_handles.append(plt.Rectangle((0,0), 1, 1, color=sex_colors[i], label=sex))
    ax.legend(handles=legend_handles, loc='upper right', bbox_to_anchor=(1, 1.2), ncol=len(sex_categories), fontsize=8)

ax_idx = 6
ax = axes[ax_idx]
marker_genes = ['SOX2',"TTR",'EOMES','GAS1','DLL1','ASCL1','DLL3','EGFR','CHGB','MKI67',
                 "SLC17A6", "SLC17A7",'UNC5D','POU3F2','CUX2','RORB','BCL11B','TBR1' , 'FOXP2','GPR39','GRIK1','ZIC4','TAFA4','NWD2',
                 'GAD1','NR2F2','SP8','LHX6','PAX7','SOX14','SIX3','DRD1','DRD2','ADARB2',
                 'GFAP','OLIG1','GPR17','MBP','TMEM119','RELN',"FOXJ1", "CLDN5",'CALD1','COL1A1']


valid_genes = [g for g in marker_genes if g in genes]
m = []
m_names = []
for gene in valid_genes:
    gene_ix = np.where(genes == gene)[0][0]
    m.append(mean_x[:, gene_ix])
    m_names.append(f"{gene}")

if m:  
    x = np.array(m)
    percentiles = np.percentile(x, 99.9, axis=1, keepdims=True)
    percentiles[percentiles == 0] = 1 
    x_norm = x / percentiles
    mask = (x == 0)
    x_norm_masked = np.ma.masked_where(mask, x_norm)
    bg = np.zeros_like(x_norm) + 0.9
    ax.imshow(bg, vmin=0, vmax=1, cmap=plt.cm.gray, aspect="auto", 
              extent=(0, len(existing_subclasses), len(m_names), 0), alpha=1)
    im = ax.imshow(x_norm_masked, cmap="inferno_r", vmin=0, vmax=1, aspect="auto", 
                   extent=(0, len(existing_subclasses), len(m_names), 0))

ax.set_yticks(np.arange(len(m_names)) + 0.5)
ax.set_yticklabels(m_names, fontsize=9)

for spine in ['top', 'right', 'bottom', 'left']:
    ax.spines[spine].set_visible(False)
    


ax_idx = 7
ax = axes[ax_idx]
im = ax.imshow(
    time_dist_ordered.T,
    vmin=0, vmax=1,
    cmap=plt.cm.bone_r,
    aspect='auto',
    interpolation="nearest",
    origin="upper",
    extent=(0, len(existing_subclasses), len(time_order), 0)
)
ax.set_yticks(np.arange(len(time_order)) + 0.5)
ax.set_yticklabels(time_order, fontsize=8)
ax.set_ylabel("time")

ax.set_xticks(np.arange(len(existing_subclasses)) + 0.5)
ax.set_xticklabels(
    existing_subclasses,
    rotation=45,
    ha="right",
    fontsize=8
)


for i, ax in enumerate(axes):
    if i < len(axes) - 1:  
        ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
        ax.spines['bottom'].set_visible(False)
    else:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        


