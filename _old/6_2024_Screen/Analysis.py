import scanpy as sc

# Read in 10X data
scData = sc.read_10x_h5("ScreenSeq_6_28_2024.h5")
scData.var_names_make_unique()

# Remove cells with high mitochondrial counts
scData.var['mt'] = scData.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(scData, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
scData = scData[scData.obs["pct_counts_mt"] < 7, :]

# Normalize and log-transform gene counts
sc.pp.normalize_total(scData, target_sum=1e4)
sc.pp.log1p(scData)

# Cell cycle regression
cell_cycle_genes = [x.strip() for x in open('regev_lab_cell_cycle_genes.txt')]
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in scData.var_names]
print("Scoring cell cycle genes...")
sc.tl.score_genes_cell_cycle(scData, s_genes=s_genes, g2m_genes=g2m_genes)

def DoClustering(data):
    print("Finding genes...")
    # Select for genes that vary across the dataset
    sc.pp.highly_variable_genes(data, min_mean=0.0125, max_mean=3, min_disp=0.5)
    data = data[:, data.var.highly_variable]

    print("Clustering...")
    # Clustering
    sc.tl.pca(data, svd_solver='arpack')
    sc.pp.neighbors(data, n_neighbors=10, n_pcs=10)
    sc.tl.leiden(data)
    sc.tl.umap(data)
    return data

scData = DoClustering(scData)

print("here")