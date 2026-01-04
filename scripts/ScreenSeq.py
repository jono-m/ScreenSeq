import typing
import pandas
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import seaborn
from sklearn.cluster import MeanShift


class ScreenSeqConfiguration:
    def __init__(self):
        self.conditionBarcodeMapping: typing.Dict[str, str] = {}
        self.normalizationBarcodes: typing.List[str] = []

    def AddConditionBarcode(self, name: str, condition: str):
        self.conditionBarcodeMapping[name] = condition

    def AddNormalizationBarcode(self, name: str):
        self.normalizationBarcodes.append(name)


def Load10X(path: str):
    # Read in 10X data
    scData = sc.read_10x_h5(path, gex_only=False)
    scData.var_names_make_unique()

    # Compute quality
    scData.var['mt'] = scData.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(scData, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    return scData


def FilterCells(scData, maxMitochondrialPercent=None,
                minTranscripts=None, maxTranscripts=None,
                minUniqueGenes=None):
    goodIndices = scData.obs["pct_counts_mt"] >= 0
    if maxMitochondrialPercent is not None:
        goodIndices &= scData.obs["pct_counts_mt"] <= maxMitochondrialPercent
    if minTranscripts is not None:
        goodIndices &= scData.obs["total_counts"] >= minTranscripts
    if maxTranscripts is not None:
        goodIndices &= scData.obs["total_counts"] <= maxTranscripts
    if minUniqueGenes is not None:
        goodIndices &= scData.obs["n_genes_by_counts"] >= minUniqueGenes
    return scData[goodIndices, :]


def FilterGenes(scData, minNumberCells=None):
    if minNumberCells is not None:
        scData = scData[:, scData.var["n_cells_by_counts"] > minNumberCells]
    return scData


def CallConditions(scData, config: ScreenSeqConfiguration):
    conditionBarcodeCounts = scData[:, list(config.conditionBarcodeMapping.keys())].to_df()
    normalizationCounts = scData[:, config.normalizationBarcodes].to_df()
    totalCounts = conditionBarcodeCounts.sum(axis=1) + normalizationCounts.sum(axis=1)

    conditionBarcodeFractions = conditionBarcodeCounts.div(totalCounts, axis=0).fillna(0)

    calls = pandas.DataFrame(np.zeros_like(conditionBarcodeFractions, dtype=bool),
                             index=conditionBarcodeFractions.index,
                             columns=conditionBarcodeFractions.columns)

    scData.uns["clusterers"] = {}
    for condition in config.conditionBarcodeMapping.keys():
        print("Clustering " + condition + "...", end="")
        amts = np.asarray(conditionBarcodeFractions[condition])
        clusterer = MeanShift(bin_seeding=True, min_bin_freq=50, n_jobs=-1).fit(amts.reshape(-1, 1))
        print(str(len(np.unique(clusterer.labels_))) + " clusters found.")
        calls[condition] = clusterer.labels_ > 0
        scData.uns["clusterers"][condition] = clusterer

    calls = calls.astype("category")

    scData.obs[[x + "_call" for x in calls.columns]] = calls
    scData.obs[[x + "_frac" for x in calls.columns]] = conditionBarcodeFractions
    scData.obs[[x + "_abs" for x in calls.columns]] = conditionBarcodeCounts

    conditions = calls.astype(str)
    conditions = (conditions.replace("False", '')
                  .replace("True", config.conditionBarcodeMapping))
    conditions = conditions.agg(lambda x: '+'.join(filter(None, x)), axis=1).astype('category')

    scData.obs["Condition"] = conditions.replace("", "Control")
    return scData


def RemoveRareConditions(scData, cutoff):
    condition, counts = np.unique(scData.obs["Condition"], return_counts=True)
    good = scData.obs["Condition"].isin(condition[counts > cutoff])
    return scData[good, :]


def Normalize(scData):
    sc.pp.normalize_total(scData)
    sc.pp.log1p(scData)
    return scData


def CellCycleScore(scData):
    cell_cycle_genes = [x.strip() for x in open('../assets/regev_lab_cell_cycle_genes.txt')]
    s_genes = cell_cycle_genes[:43]
    g2m_genes = cell_cycle_genes[43:]
    sc.tl.score_genes_cell_cycle(scData, s_genes=s_genes, g2m_genes=g2m_genes)
    return scData


def RegressCellCycle(scData):
    sc.pp.regress_out(scData, ['S_score', 'G2M_score'])
    return scData


def Cluster(scData):
    print("Finding variable genes...")
    sc.pp.highly_variable_genes(scData)
    scData = scData[:, scData.var.highly_variable]
    print("Clustering...")
    sc.tl.pca(scData, svd_solver='arpack')
    sc.pp.neighbors(scData, n_neighbors=10, n_pcs=10)
    sc.tl.leiden(scData)
    sc.tl.umap(scData)
    return scData


def ClusterScreenSeq(scData):
    scData2 = scData[:, scData.var_names.str.startswith("ScreenSeq")]
    sc.tl.pca(scData2, svd_solver='arpack')
    sc.pp.neighbors(scData2, n_neighbors=10, n_pcs=10, key_added="screenseq")
    sc.tl.umap(scData, neighbors_key="screenseq")



def CompareGene(scData, gene, cytokine, config: ScreenSeqConfiguration):
    barcode = [x for x, y in config.conditionBarcodeMapping.items() if y == cytokine][0]
    barcode = barcode + "_call"
    geneData = scData[:, gene].to_df().join(scData.obs[barcode])
    geneData[barcode] = np.where(geneData[barcode] == True, "+" + cytokine, "Control")
    seaborn.violinplot(geneData, x=gene, y=barcode, fill=False)
    plt.xlabel(gene + " expression")
    plt.ylabel("")
    plt.show()
