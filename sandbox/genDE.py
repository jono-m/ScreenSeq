import pandas as pd
from pydeseq2.dds import DeseqDataSet, DefaultInference

countPath = r"../data/Sequencing/BulkSeq_2025_03_06/counts.csv"
metadataPath = r"../data/Sequencing/BulkSeq_2025_03_06/metadata.csv"
counts = pd.read_csv(countPath, index_col=0)
metadata = pd.read_csv(metadataPath, index_col=0)

ligand = "TLR-1"
samplesToKeep = metadata["Ligand"].isin([ligand, "Control"]) & metadata["Dose"].isin(["high", "none"])

metadata = metadata[samplesToKeep]
counts = counts[samplesToKeep]

genes_to_keep = counts.columns[counts.sum(axis=0) >= 10]
counts = counts[genes_to_keep]

inference = DefaultInference(n_cpus=None)
dds = DeseqDataSet(counts=counts,metadata=metadata,
                   design="~Ligand", refit_cooks=True,inference=inference)
dds.deseq2()

import pickle
with open("./DE.pkl", "wb+") as f:
    pickle.dump(dds, f)