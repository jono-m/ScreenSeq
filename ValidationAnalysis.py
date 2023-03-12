import pandas
import collections
import tables
import scipy.sparse
import seaborn
import matplotlib.pyplot as plt
import numpy as np


def get_matrix_from_h5(filename):
    CountMatrix = collections.namedtuple('CountMatrix', ['feature_ref', 'barcodes', 'matrix'])
    with tables.open_file(filename, 'r') as f:
        mat_group = f.get_node(f.root, 'matrix')
        barcodes = f.get_node(mat_group, 'barcodes').read()
        data = getattr(mat_group, 'data').read()
        indices = getattr(mat_group, 'indices').read()
        indptr = getattr(mat_group, 'indptr').read()
        shape = getattr(mat_group, 'shape').read()
        matrix = scipy.sparse.csc_matrix((data, indices, indptr), shape=shape)

        feature_ref = {}
        feature_group = f.get_node(mat_group, 'features')
        feature_ids = getattr(feature_group, 'id').read()
        feature_names = getattr(feature_group, 'name').read()
        feature_types = getattr(feature_group, 'feature_type').read()
        feature_ref['id'] = feature_ids
        feature_ref['name'] = feature_names
        feature_ref['feature_type'] = feature_types
        tag_keys = getattr(feature_group, '_all_tag_keys').read()
        for key in tag_keys:
            key = key.decode("utf-8")
            feature_ref[key] = getattr(feature_group, key).read()

        return CountMatrix(feature_ref, barcodes, matrix)


def mapBCName(a):
    d = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return a[:7] + d[a[7]] + d[a[8]] + a[9:]


def Analyze():
    ssData = get_matrix_from_h5(r"Data\3_8_2023_Validation\filtered_feature_bc_matrix.h5")
    ssBarcodeIndices = np.where(np.char.startswith(ssData[0]["name"].astype(str), "ScreenSeq"))[0]
    ssBarcodes = list(ssData[0]['name'][ssBarcodeIndices].astype(str))
    ssBarcodes = ["SS" + x[9:] for x in ssBarcodes]
    cellBarcodes = np.asarray([x[:-2] for x in ssData[1].astype(str)])
    ssData = ssData[2][ssBarcodeIndices, :].toarray().transpose()
    ssData = pandas.DataFrame(data=np.concatenate([cellBarcodes[:, None], ssData], axis=1),
                              columns=["CellBC"] + ssBarcodes).set_index("CellBC").astype(int)

    # ssData = pandas.read_csv(r"Data\3_8_2023_Validation\ssCounts_HD2.csv")
    # ssData["CellBC"] = ssData["CellBC"].map(mapBCName)
    # ssData = ssData.set_index("CellBC")
    # ssBarcodes = ssData.columns

    ssData["SSTot"] = ssData.sum(axis=1)
    ssData[[bc + "_norm" for bc in ssBarcodes]] = ssData[ssBarcodes].div(ssData["SSTot"], axis=0)

    cellData = get_matrix_from_h5(r"Data\3_8_2023_Validation\filtered_feature_bc_matrix_cDNA.h5")
    mtGeneNumbers = np.where(np.char.startswith(cellData[0]["name"].astype(str), "MT"))[0]
    mtGeneCounts = cellData[2][mtGeneNumbers, :].toarray().sum(axis=0)
    totalUMIs = np.asarray(cellData[2].sum(axis=0))[0]
    uniqueGene = np.asarray((cellData[2] > 0).sum(axis=0))[0]
    cellBarcodes = np.asarray([x[:-2] for x in cellData[1].astype(str)])
    cellData = np.asarray([cellBarcodes, mtGeneCounts, totalUMIs, uniqueGene])
    cellColumns = ["CellBC", "Mitochondrial UMIs", "Total UMIs", "Unique Features"]
    cellData = pandas.DataFrame(data=cellData.transpose(), columns=cellColumns).set_index("CellBC").astype(int)
    cellData["Percent MT"] = cellData["Mitochondrial UMIs"] / cellData["Total UMIs"]

    fullData = pandas.concat([ssData, cellData], axis=1)
    liveCells = fullData[(fullData["Unique Features"] > 4000) &
                         (fullData["Percent MT"] < 0.1) &
                         (fullData["Total UMIs"] > 18000) &
                         (fullData["Total UMIs"] < 80000)]
    fullPass = liveCells
    fig = plt.figure()
    # plt.subplot(2, 2, 1)
    # seaborn.kdeplot(data=fullPass, x="SS4_norm", y="SS5_norm", color="black")
    # seaborn.scatterplot(data=fullPass, x="SS4_norm", y="SS5_norm", color="black")
    #
    # plt.subplot(2, 2, 2)
    # seaborn.kdeplot(data=fullPass, x="SS6_norm", y="SS7_norm", color="orange")
    # seaborn.scatterplot(data=fullPass, x="SS6_norm", y="SS7_norm", color="orange")

    # ax = fig.add_subplot(2, 2, 3, projection="3d")
    ax = fig.add_subplot(projection="3d")
    ax.scatter(xs=fullPass["SS8_norm"], ys=fullPass["SS9_norm"], zs=fullPass["SS10_norm"], color='red')
    ax.set_xlabel("SS8")
    ax.set_ylabel("SS9")
    ax.set_zlabel("SS10")

    # plt.subplot(2, 2, 4)
    # seaborn.kdeplot(data=fullPass, x="SS11_norm", y="SS10_norm", color="green")
    # seaborn.scatterplot(data=fullPass, x="SS11_norm", y="SS10_norm", color="green")
    plt.show()


Analyze()
