import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import harmonypy as hm


def run_preprocess(
    adata,
    source_layer="counts",
    target_sum=1e4,
    scale_max=10,
    n_comps=20,
    svd_solver="arpack",
):
    """Set counts to X, normalize/log, scale, and run PCA."""
    adata.X = adata.layers[source_layer].copy()
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, max_value=scale_max)
    sc.tl.pca(adata, svd_solver=svd_solver, n_comps=n_comps)
    return adata


def run_harmony(
    adata,
    batch_key,
    input_pca_key="X_pca",
    output_key="X_pca_harmony",
    max_iter_harmony=30,
    theta=1,
    verbose=True,
):
    """Run Harmony from an existing PCA embedding and store corrected PCs."""
    data_mat = adata.obsm[input_pca_key]
    meta_data = adata.obs
    vars_use = [batch_key]

    print(f"Starting Harmony on {data_mat.shape} matrix...")
    ho = hm.run_harmony(
        data_mat,
        meta_data,
        vars_use,
        max_iter_harmony=max_iter_harmony,
        verbose=verbose,
        theta=theta,
    )

    res = ho.Z_corr
    if hasattr(res, "cpu"):
        res = res.cpu().numpy()
    elif hasattr(res, "numpy"):
        res = res.numpy()

    if res.ndim != 2:
        raise ValueError(f"Harmony output is not 2D: shape={res.shape}")

    if res.shape[1] == adata.shape[0]:
        res = res.T

    print(f"Harmony Output Shape: {res.shape}")
    if res.shape[0] != adata.shape[0]:
        raise ValueError(f"Dimension mismatch. AnnData: {adata.shape[0]}, Harmony: {res.shape[0]}")

    adata.obsm[output_key] = res
    print(f"Integration successful. Added '{output_key}' to obsm.")
    return adata


def run_umapclust(
    adata,
    use_rep="X_pca_harmony",
    n_neighbors=15,
    n_pcs=20,
    leiden_resolution=0.1,
    leiden_key="leiden",
):
    """Build graph, compute UMAP, and cluster with Leiden."""
    sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=leiden_resolution, key_added=leiden_key)
    return adata


def run_deg_prep(adata, source_layer="counts", target_sum=1e4):
    """Set counts to X and apply log-normalization for DE testing."""
    adata.X = adata.layers[source_layer].copy()
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    return adata


def run_deg(adata, groupby="leiden", method="wilcoxon", use_raw=False, pts=True):
    """Run Scanpy marker ranking and store results in adata.uns."""
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        method=method,
        use_raw=use_raw,
        pts=pts,
    )
    return adata


def build_top_deg_table(
    adata,
    top_n=25,
    pval_adj_max=0.05,
    min_log2fc=None,
    result_key="rank_genes_groups",
):
    """
    Build a per-cluster DEG summary table from adata.uns['rank_genes_groups'].

    Returns columns: cluster, gene, log2fc, pvals_adj, pct_in_group.
    """
    if result_key not in adata.uns:
        raise KeyError(f"'{result_key}' not found in adata.uns. Run run_deg(...) first.")

    result = adata.uns[result_key]
    groups = result["names"].dtype.names
    rows = []

    for group in groups:
        df = pd.DataFrame(
            {
                "gene": result["names"][group],
                "log2fc": result["logfoldchanges"][group],
                "pvals_adj": result["pvals_adj"][group],
                "pct_in_group": result["pts"][group].values if "pts" in result else np.nan,
            }
        )

        filt = df["pvals_adj"] < pval_adj_max
        if min_log2fc is not None:
            filt &= df["log2fc"] >= min_log2fc

        df = df.loc[filt].sort_values("log2fc", ascending=False).head(top_n).copy()
        df["cluster"] = group
        rows.append(df)

    if not rows:
        return pd.DataFrame(columns=["cluster", "gene", "log2fc", "pvals_adj", "pct_in_group"])

    out = pd.concat(rows, ignore_index=True)
    return out[["cluster", "gene", "log2fc", "pvals_adj", "pct_in_group"]]


def map_clusters(
    adata,
    cluster_map,
    source_cluster_key="leiden",
    target_celltype_key="cell_type",
    strict=True,
):
    """
    Map cluster IDs to labels.

    cluster_map should be a dict with cluster IDs as strings, e.g.:
    {'0': 'Macrophages', '1': 'T cells', '2': 'B cells'}
    """
    adata.obs[target_celltype_key] = adata.obs[source_cluster_key].astype(str).map(cluster_map)
    if strict and adata.obs[target_celltype_key].isnull().any():
        missing = adata.obs.loc[adata.obs[target_celltype_key].isnull(), source_cluster_key].astype(str).unique()
        raise ValueError(f"Unmapped clusters in {source_cluster_key}: {missing}")
    return adata


def sanitize_for_h5ad(adata):
    """Replace '/' with '_' in obs columns and uns keys for HDF5 compatibility."""
    adata.obs.columns = [c.replace("/", "_") for c in adata.obs.columns]
    for key in list(adata.uns.keys()):
        if "/" in key:
            adata.uns[key.replace("/", "_")] = adata.uns.pop(key)
    return adata


def calculate_density_grid(
    adata_obj,
    compartments,
    resolution=50,
    celltype_col="celltype_imm_agg",
):
    """
    Calculate cell density per compartment using a spatial grid approximation.

    compartments: dict mapping displayed compartment name -> boolean column in .obs
    example: {'Glomerulus': 'is_in_glom', 'Periglomerular': 'is_in_periglom'}
    """
    results = []
    for slide_id in adata_obj.obs["slide"].unique():
        subset = adata_obj.obs[adata_obj.obs["slide"] == slide_id].copy()

        x_min, x_max = subset["x_centroid"].min(), subset["x_centroid"].max()
        y_min, y_max = subset["y_centroid"].min(), subset["y_centroid"].max()
        x_bins = np.arange(x_min, x_max + resolution, resolution)
        y_bins = np.arange(y_min, y_max + resolution, resolution)

        subset["x_bin"] = np.digitize(subset["x_centroid"], x_bins)
        subset["y_bin"] = np.digitize(subset["y_centroid"], y_bins)
        subset["bin_id"] = subset["x_bin"].astype(str) + "_" + subset["y_bin"].astype(str)

        slide_areas = {}
        for name, col in compartments.items():
            comp_cells = subset[subset[col].isin([True, "True", 1, "1"])]
            n_occupied_bins = comp_cells["bin_id"].nunique()
            area_um2 = n_occupied_bins * (resolution ** 2)
            area_mm2 = area_um2 / 1_000_000
            slide_areas[name] = area_mm2
            print(f"Slide {slide_id} - {name} Area: {area_mm2:.4f} mm^2")

        for name, col in compartments.items():
            comp_cells = subset[subset[col].isin([True, "True", 1, "1"])]
            counts = comp_cells[celltype_col].value_counts()
            denom_area = slide_areas[name]

            if denom_area > 0:
                densities = counts / denom_area
                for cell_type, dens in densities.items():
                    results.append(
                        {
                            "Slide": slide_id,
                            "Compartment": name,
                            "Cell Type": cell_type,
                            "Count": counts[cell_type],
                            "Area_mm2": denom_area,
                            "Density_cells_mm2": dens,
                        }
                    )

    return pd.DataFrame(results)


def plot_density(
    adata_obj,
    group_by,
    split_by,
    replicate_by=None,
    normalize="index",
    figsize=(9, 4),
    dpi=200,
    fliersize=0.5,
    rotation=45,
):
    """
    Plot .obs composition as boxplot using group/split variables.

    If replicate_by is provided, crosstab is computed per replicate and boxplots
    show distribution across replicates.
    """
    needed = [group_by, split_by] + ([replicate_by] if replicate_by else [])
    for col in needed:
        if col not in adata_obj.obs.columns:
            raise KeyError(f"Column '{col}' not found in adata.obs")

    obs_df = adata_obj.obs[needed].dropna().copy()

    if replicate_by:
        frames = []
        for rep, sub in obs_df.groupby(replicate_by):
            ctab = pd.crosstab(sub[group_by], sub[split_by], normalize=normalize)
            long_df = ctab.reset_index().melt(id_vars=group_by, var_name=split_by, value_name="value")
            long_df[replicate_by] = rep
            frames.append(long_df)
        plot_df = pd.concat(frames, ignore_index=True)
    else:
        ctab = pd.crosstab(obs_df[group_by], obs_df[split_by], normalize=normalize)
        plot_df = ctab.reset_index().melt(id_vars=group_by, var_name=split_by, value_name="value")

    with plt.rc_context({"figure.figsize": figsize, "figure.dpi": dpi}):
        ax = sns.boxplot(plot_df, x=group_by, y="value", hue=split_by, fliersize=fliersize)
        plt.setp(ax.get_xticklabels(), rotation=rotation)
        plt.tight_layout()
        plt.show()
    return plot_df


def get_bin(coords, min_val, res):
    """Convert coordinates to integer grid-bin index."""
    return ((coords - min_val) // res).astype(int)
