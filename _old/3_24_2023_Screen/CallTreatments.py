import pandas
import tables
import scipy.sparse
import seaborn
import matplotlib.pyplot as plt
import numpy as np
from sklearn.cluster import MeanShift


def loadH5Data(filename):
    with tables.open_file(filename, 'r') as f:
        mat_group = f.get_node(f.root, 'matrix')
        cellBarcodes = np.asarray(
            [x[:-2] for x in f.get_node(mat_group, 'barcodes').read().astype(str)])
        data = getattr(mat_group, 'data').read()
        indices = getattr(mat_group, 'indices').read()
        indptr = getattr(mat_group, 'indptr').read()
        shape = getattr(mat_group, 'shape').read()
        matrix = scipy.sparse.csc_matrix((data, indices, indptr), shape=shape)

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


def ComputeCellQualityMetrics(features, cellBarcodes, expressionMatrix):
    ssBarcodeIndices = np.char.startswith(features["name"], "ScreenSeq")
    mitochondialIndices = np.char.startswith(features["name"], "MT")
    mitochondrialUMIs = expressionMatrix[:, mitochondialIndices].sum(axis=1)
    totalUMIs = expressionMatrix[:, ~ssBarcodeIndices].sum(axis=1)
    uniqueGenes = (expressionMatrix[:, ~ssBarcodeIndices] > 0).sum(axis=1)
    return mitochondrialUMIs, uniqueGenes, totalUMIs


def RemoveBadCells(features, cellBarcodes, expressionMatrix):
    _, uniqueGenes, totalUMIs = ComputeCellQualityMetrics(features, cellBarcodes, expressionMatrix)
    goodCells = (uniqueGenes > 4000) & (totalUMIs < 40000) & (totalUMIs > 10000)
    return features, cellBarcodes[goodCells], expressionMatrix[goodCells]


def LoadSSData():
    # Load 10X feature expression matrix
    features, cellBarcodes, expressionMatrix = RemoveBadCells(
        *loadH5Data(r"raw_feature_bc_matrix.h5"))

    # Grab Screen-seq barcode UMI counts
    ssBarcodeIndices = np.char.startswith(features["name"], "ScreenSeq")
    ssBarcodeNames = list(features['name'][ssBarcodeIndices])
    ssData = expressionMatrix[:, ssBarcodeIndices]

    # Only use barcodes 4-11 and N.
    indicesToUse = [3, 4, 5, 6, 7, 8, 9, 10, 16]
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
    ssData = ssData[:, :-1]

    # Concatenate to a Pandas frame. Assign the drug names.
    ssData = pandas.DataFrame(data=ssData,
                              columns=["Tanespimycin", "Cytarabine", "Dasatinib",
                                       "Omacetaxine", "Hydroxyurea",
                                       "Imatinib", "Ruxolitinib", "Vorinostat"],
                              index=cellBarcodes).astype(float)
    ssData = ssData[sorted(ssData.columns)] * 100

    return ssData


def NormalizeKDEPlotHeight(ax):
    for collection in ax.collections:
        yCoords = collection.get_paths()[0].vertices[:, 1]
        yCoords = (yCoords - np.min(yCoords)) / (np.max(yCoords) - np.min(yCoords))
        collection.get_paths()[0].vertices[:, 1] = yCoords
    for line in ax.lines:
        yCoords = line.get_ydata()
        yCoords = (yCoords - np.min(yCoords)) / (np.max(yCoords) - np.min(yCoords))
        line.set_ydata(yCoords)


def CallSSData(data):
    calls = []
    clusterers = []
    for drug in data.columns:
        barcodeAmounts = data[drug]

        clusterer = MeanShift(bin_seeding=True, min_bin_freq=round(len(barcodeAmounts) / 10)).fit(
            np.asarray(barcodeAmounts).reshape(-1, 1))
        calls.append((drug, clusterer.labels_ > 0))
        clusterers.append(clusterer)

    df = pandas.DataFrame(data=np.asarray([response for drug, response in calls]).T,
                          index=data.index,
                          columns=[drug for drug, response in calls])
    return df, clusterers


def PlotCalls(data, calls, clusterers):
    fig, axs = plt.subplots(2, 4, sharex="col", sharey="row", gridspec_kw=dict(wspace=0, hspace=0))

    for drug, ax, color, clusterer in zip(data.columns, axs.flatten(),
                                          seaborn.color_palette("husl", n_colors=len(data.columns)),
                                          clusterers):
        plt.sca(ax)

        barcodeAmounts = data[drug]
        kdeAxis = seaborn.kdeplot(barcodeAmounts, ax=ax, fill=False, color=[0, 0, 0, 0.2],
                                  alpha=0.2,
                                  clip=[0, 100])
        NormalizeKDEPlotHeight(kdeAxis)

        vertices = np.stack(kdeAxis.lines[0].get_data(), axis=-1)
        verticesClustered = clusterer.predict(vertices[:, 0:1])
        splitCoords = [0] + list(np.where(np.diff(verticesClustered) != 0)[0] + 1) + [-1]
        segments = zip(splitCoords[0:], splitCoords[1:])
        for start, end in segments:
            cluster = verticesClustered[start]
            end = min(-1, end + 1)
            pts = vertices[start:end, :]
            ax.fill_between(x=pts[:, 0], y1=pts[:, 1],
                            color=[0.7, 0.7, 0.7] if cluster == 0 else color,
                            alpha=0.5)

        plt.scatter(x=data.loc[~calls[drug], drug],
                    y=np.random.uniform(-1, 1, np.count_nonzero(~calls[drug])) * 0.025 - 0.05,
                    color=[0.7, 0.7, 0.7], s=0.25)
        plt.scatter(x=data.loc[calls[drug], drug],
                    y=np.random.uniform(-1, 1, np.count_nonzero(calls[drug])) * 0.025 - 0.05,
                    color=color, s=0.25)
        ax.set_yticks([])
        ax.set_ylabel(None)
        ax.set_xlabel(None)
        ax.set_ylim([-0.1, 1.1])
        ax.set_xlim([0, 40])
        ax.set_xticks([0, 10, 20, 30])
        plt.text(np.mean(ax.get_xlim()), 0.9, drug, ha="center", color=color)

    # Turn off axis lines and ticks of the big subplot
    fig.supxlabel("Drug barcode relative abundance (%)")
    fig.supylabel("Normalized cell counts")
    plt.tight_layout()
    plt.show()


def Run():
    data = LoadSSData()
    calledData, clusterers = CallSSData(data)
    calledData.to_csv("Calls_raw.csv")
    PlotCalls(data, calledData, clusterers)


Run()
