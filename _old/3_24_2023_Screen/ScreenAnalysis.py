import pandas
import numpy as np
import tables
import scipy
import seaborn
import matplotlib.pyplot as plt


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


def ComputeCellQualityMetrics(features, expressionMatrix):
    ssBarcodeIndices = np.char.startswith(features["name"], "ScreenSeq")
    mitochondialIndices = np.char.startswith(features["name"], "MT")
    mitochondrialUMIs = expressionMatrix[:, mitochondialIndices].sum(axis=1)
    totalUMIs = expressionMatrix[:, ~ssBarcodeIndices].sum(axis=1)
    uniqueGenes = (expressionMatrix[:, ~ssBarcodeIndices] > 0).sum(axis=1)
    return mitochondrialUMIs, uniqueGenes, totalUMIs


def RemoveBadCells(features, cellBarcodes, expressionMatrix):
    _, uniqueGenes, totalUMIs = ComputeCellQualityMetrics(features, expressionMatrix)
    goodCells = (uniqueGenes > 4000) & (totalUMIs < 40000) & (totalUMIs > 10000)
    return features, cellBarcodes[goodCells], expressionMatrix[goodCells]


def Run():
    treatmentData = pandas.read_csv("Calls.csv", index_col=0)
    features, cellBarcodes, expressionMatrix = loadH5Data("ScreenSeq_6_28_2024.h5")
    print("Loaded treatment data for %d cells." % len(treatmentData))

    metrics = ComputeCellQualityMetrics(features, expressionMatrix)
    df = pandas.DataFrame(data=np.asarray(metrics).T,
                          columns=["Mitochondrial UMIs", "Unique Genes",
                                   "Total UMIs"],
                          index=cellBarcodes)
    df.index.rename("Cell", inplace=True)

    data = treatmentData.join(df)
    data["MitochondrialFrac"] = data["Mitochondrial UMIs"] / data["Total UMIs"]

    data = data.sort_values(by="MitochondrialFrac")
    seaborn.heatmap(data=data[list(treatmentData.columns) + ["MitochondrialFrac"]])
    plt.show()

Run()
