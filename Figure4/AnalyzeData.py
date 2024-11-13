import scanpy as sc

# Read in 10X data
scData = sc.read_10x_h5("../data/Sequencing/screenSeq_9_2024.h5", gex_only=False)
scData.var_names_make_unique()

# Compute quality
scData.var['mt'] = scData.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(scData, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)