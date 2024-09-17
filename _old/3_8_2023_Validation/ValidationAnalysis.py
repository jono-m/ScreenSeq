import pandas
import collections
import tables
import scipy.sparse
import seaborn
import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import linkage as Linkage


def loadH5Data(filename):
    with tables.open_file(filename, 'r') as f:
        mat_group = f.get_node(f.root, 'matrix')
        cellBarcodes = np.asarray(
            [x[:-2] for x in f.get_node(mat_group, 'barcodes').read().astype(str)])
        data = getattr(mat_group, 'data').read()
        indices = getattr(mat_group, 'indices').read()
        indptr = getattr(mat_group, 'indptr').read()
        shape = getattr(mat_group, 'shape').read()
        matrix = scipy.sparse.csc_matrix((data, indices, indptr), shape=shape).toarray().transpose()

        features = {}
        feature_group = f.get_node(mat_group, 'features')
        feature_ids = getattr(feature_group, 'id').read().astype(str)
        feature_names = getattr(feature_group, 'name').read().astype(str)
        feature_types = getattr(feature_group, 'feature_type').read().astype(str)
        features['id'] = feature_ids
        features['name'] = feature_names
        features['feature_type'] = feature_types
        tag_keys = getattr(feature_group, '_all_tag_keys').read()
        for key in tag_keys:
            key = key.decode("utf-8")
            features[key] = getattr(feature_group, key).read().astype(str)

        return features, cellBarcodes, matrix


def RemoveBadCells(features, cellBarcodes, expressionMatrix):
    ssBarcodeIndices = np.char.startswith(features["name"], "ScreenSeq")
    mitochondialIndices = np.char.startswith(features["name"], "MT")
    mitochondrialUMIs = expressionMatrix[:, mitochondialIndices].sum(axis=1)
    totalUMIs = expressionMatrix[:, ~ssBarcodeIndices].sum(axis=1)
    uniqueGenes = (expressionMatrix[:, ~ssBarcodeIndices] > 0).sum(axis=1)

    goodCells = (uniqueGenes > 4000) & (mitochondrialUMIs/totalUMIs < 0.1) & (totalUMIs < 60000) & (totalUMIs > 20000)
    return features, cellBarcodes[goodCells], expressionMatrix[goodCells]


def LoadSSData():
    # Load 10X feature expression matrix
    features, cellBarcodes, expressionMatrix = RemoveBadCells(
        *loadH5Data(r"filtered_feature_bc_matrix_REDO.h5"))

    # Grab Screen-seq barcode UMI counts
    ssBarcodeIndices = np.char.startswith(features["name"], "ScreenSeq")
    ssBarcodeNames = features['name'][ssBarcodeIndices]
    ssData = expressionMatrix[:, ssBarcodeIndices]

    # Only use barcodes 4-11 and N.
    indicesToUse = [3, 4, 5, 6, 7, 8, 9, 10, 16]
    ssBarcodeNames = list(ssBarcodeNames[indicesToUse])
    ssData = ssData[:, indicesToUse]

    # Drop cells with abnormally high or low numbers of total barcodes
    totalBarcodes = ssData.sum(axis=1)
    zScore = np.abs(totalBarcodes - totalBarcodes.mean()) / totalBarcodes.std()
    goodCells = (zScore <= 3) & (totalBarcodes > 0)

    ssData = ssData[goodCells, :]
    cellBarcodes = cellBarcodes[goodCells, None]

    # Normalize per-cell.
    ssData = ssData / ssData.sum(axis=1, keepdims=True)

    # Drop the normalization barcode
    ssBarcodeNames = ssBarcodeNames[:-1]
    ssData = ssData[:, :-1]

    # Concatenate to a Pandas frame.
    ssData = pandas.DataFrame(np.concatenate([cellBarcodes, ssData], axis=1),
                              columns=["Cell"] + ssBarcodeNames)
    ssData = ssData.set_index("Cell").astype(float)
    ssData = ssData * 100

    return ssData


def Analyze():
    ssData = LoadSSData()
    ssData = ssData.loc[:, [c for c in ssData.columns if c != "ScreenSeq11"]]

    col_linkage = Linkage(np.asarray(ssData).T, method="average", optimal_ordering=True)
    g = seaborn.clustermap(ssData, standard_scale=1,
                           xticklabels=["MonoA", "MonoB", "Duo-1", "Duo-2", "Trio-2",
                                        "Trio-1", "Trio-3"],
                           yticklabels=[],
                           dendrogram_ratio=(0.1, 0.1),
                           figsize=(8, 7),
                           cbar_pos=(1, 0.2, 0.1, 0.6),
                           col_linkage=col_linkage)
    g.fig.subplots_adjust(right=0.8, bottom=0.1)
    g.ax_cbar.set_position([0.85, 0.1, 0.02, 0.78])
    g.ax_cbar.set_ylabel("Normalized barcode abundance", labelpad=10)
    g.ax_heatmap.set_xlabel("Barcode")
    g.ax_heatmap.set_ylabel("Cell (n=%d)" % len(ssData))
    plt.show()


Analyze()
