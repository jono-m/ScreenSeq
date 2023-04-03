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


def LoadSSData():
    # Load 10X feature expression matrix
    features, cellBarcodes, expressionMatrix = loadH5Data(r"filtered_feature_bc_matrix.h5")

    # Grab Screen-seq barcode UMI counts
    ssBarcodeIndices = np.char.startswith(features["name"], "ScreenSeq")
    ssBarcodeNames = list(features['name'][ssBarcodeIndices])
    ssData = expressionMatrix[:, ssBarcodeIndices]

    # Only use barcodes 4-11 and N.
    indicesToUse = [3, 4, 5, 6, 7, 8, 9, 10, 16]
    ssBarcodeNames = [ssBarcodeNames[indexToUse] for indexToUse in indicesToUse]
    ssData = ssData[:, indicesToUse]

    # Drop cells that have no barcodes detected
    goodCells = ssData.sum(axis=1) > 0
    ssData = ssData[goodCells, :]
    cellBarcodes = cellBarcodes[goodCells, None]

    # Normalize per-cell.
    ssData = ssData / ssData.sum(axis=1, keepdims=True)

    # Drop the normalization barcode
    ssData = ssData[:, :-1]

    # Concatenate to a Pandas frame. Assign the drug names.
    ssData = pandas.DataFrame(np.concatenate([cellBarcodes, ssData], axis=1),
                              columns=["Cell", "Tanespimycin", "Cytarabine", "Dasatinib",
                                       "Homoharringtonine", "Hydroxyurea",
                                       "Imatinib", "Ruxolitinib", "Vorinostat"])
    ssData = ssData.set_index("Cell").astype(float)
    ssData = ssData.loc[(ssData <= 0.4).all(axis=1)]

    return ssData


# Cluster cells to different treatment groups based on barcode abundances
def ClusterTreatments(ssData):
    numDrugs = len(ssData.columns)
    centroidPermutations = np.asarray([[int(x) for x in format(i, "0%db" % numDrugs)] for i in range(2 ** numDrugs)])
    centroidLo = np.asarray(ssData.quantile(0.25))
    centroidHi = np.asarray(ssData.quantile(0.75))
    centroids = (centroidPermutations * (centroidHi-centroidLo))+centroidLo

    distanceMatrix = scipy.spatial.distance_matrix(np.asarray(ssData), centroids)
    closestCentroidIndices = np.argmin(distanceMatrix, axis=1)
    closestCentroids = centroids[closestCentroidIndices]
    distanceToClosest = distanceMatrix.min(axis=1, keepdims=True)
    return pandas.DataFrame(index=ssData.index, data=np.concatenate([closestCentroidIndices[:, None], closestCentroids, distanceToClosest], axis=1),
                            columns=["Treatment ID"] + list(ssData.columns) + ["Distance"])


data = LoadSSData()
clustered = ClusterTreatments(data)
data["TreatmentID"] = clustered["Treatment ID"]
# data = data.loc[clustered["Distance"] <= 0.1]

data = data.sort_values("TreatmentID", ascending=False)
g = seaborn.clustermap(data.iloc[:, :-1], standard_scale=1,
                       yticklabels=[],
                       figsize=(8, 7),
                       cbar_pos=(1, 0.2, 0.1, 0.6),
                       col_cluster=False,
                       row_cluster=False,
                       cmap="viridis")
g.fig.subplots_adjust(right=0.8, bottom=0.1)
g.ax_cbar.set_position([0.85, 0.1, 0.02, 0.78])
g.ax_cbar.set_ylabel("Normalized barcode abundance", labelpad=10)
g.ax_heatmap.tick_params(axis='x', rotation=30)
g.ax_heatmap.set_xlabel("Barcode")
g.ax_heatmap.set_ylabel("Cell (n=%d)" % len(data))
plt.show()
