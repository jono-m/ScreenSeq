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


def Analyze():
    # Load 10X feature expression matrix
    features, cellBarcodes, expressionMatrix = loadH5Data(r"filtered_feature_bc_matrix_REDO.h5")

    # Extract Screen-seq barcodes.
    ssBarcodeIndices = np.char.startswith(features["name"], "ScreenSeq")
    ssBarcodeNames = list(features['name'][ssBarcodeIndices])
    cellExpressionData = expressionMatrix[:, ssBarcodeIndices]
    cellExpressionData = np.concatenate([cellBarcodes[:, None], cellExpressionData], axis=1)
    cellExpressionData = pandas.DataFrame(cellExpressionData, columns=["Cell"] + ssBarcodeNames)
    cellExpressionData = cellExpressionData.set_index("Cell").astype(int)

    # Normalize Screen-seq barcodes.
    cellExpressionData["ScreenSeqTotal"] = cellExpressionData.sum(axis=1)
    cellExpressionData[[bc + "_norm" for bc in ssBarcodeNames]] = cellExpressionData[
        ssBarcodeNames].divide(cellExpressionData["ScreenSeqTotal"], axis=0)

    # Count the number of mitochondrial UMIs per cell, the number of unique genes per cell, and the
    # total UMI counts per cell.
    mitochondrialIndices = np.char.startswith(features["name"], "MT")
    cellExpressionData["Mitochondrial UMIs"] = \
        expressionMatrix[:, mitochondrialIndices].sum(axis=1, keepdims=True)

    geneIndices = ~(np.char.startswith(features["name"], "ScreenSeq") |
                    np.char.startswith(features["name"], "MT"))
    geneExpressionMatrix = expressionMatrix[:, geneIndices]
    cellExpressionData["Gene UMIs"] = geneExpressionMatrix.sum(axis=1, keepdims=True)
    cellExpressionData["Genes"] = (geneExpressionMatrix > 0).sum(axis=1)

    # Ignore cells with low screen seq numbers, genes, UMIs, or high mitochondrial UMIs.
    # cellExpressionData = cellExpressionData.loc[
    #     (cellExpressionData["Genes"] > 2000) &
    #     (cellExpressionData["Gene UMIs"] > 12000) &
    #     (cellExpressionData["ScreenSeqTotal"] > 0) &
    #     (cellExpressionData["Mitochondrial UMIs"] < 3000) &
    #     (cellExpressionData["ScreenSeqN_norm"] < 0.95)
    # ]
    cellExpressionData = cellExpressionData.loc[
        (cellExpressionData["Genes"] > 2000) &
        (cellExpressionData["Gene UMIs"] > 12000) &
        (cellExpressionData["ScreenSeqTotal"] > 100) &
        (cellExpressionData["Mitochondrial UMIs"] < 3000) &
        (cellExpressionData["ScreenSeqN_norm"] < 0.93)
        ]
    ssData = cellExpressionData[["%s_norm" % s for s in ssBarcodeNames[3:-7]]].copy()

    l2 = Linkage(np.asarray(ssData).T, method="average", optimal_ordering=True)
    g = seaborn.clustermap(ssData, standard_scale=1,
                           xticklabels=["MonoA", "MonoB", "Duo-1", "Duo-2", "Trio-2",
                                        "Trio-1", "Trio-3"],
                           yticklabels=[],
                           col_colors=["#AA2222", "#22AA22", "#2222AA", "#2222AA", "#AA22AA", "#AA22AA",
                            "#AA22AA"],
                           dendrogram_ratio=(0.1, 0.1),
                           figsize=(8, 7),
                           cbar_pos=(1, 0.2, 0.1, 0.6),
                           col_linkage=l2)
    g.fig.subplots_adjust(right=0.8, bottom=0.1)
    g.ax_cbar.set_position([0.85, 0.1, 0.02, 0.78])
    g.ax_cbar.set_ylabel("Normalized barcode abundance", labelpad=10)
    g.ax_heatmap.set_xlabel("Barcode")
    g.ax_heatmap.set_ylabel("Cell (n=%d)" % len(ssData))
    plt.show()


Analyze()
