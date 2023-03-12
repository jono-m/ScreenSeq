import pandas
import collections
import tables
import scipy.sparse
import seaborn
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

def Analyze():
    ssData = pandas.read_csv(r"Data\3_8_2023_Validation\ssCounts_COMP3_12.csv")
    cellData = get_matrix_from_h5(r"Data\3_8_2023_Validation\filtered_feature_bc_matrix.h5")
    mtGeneNumbers = np.where(np.char.startswith(cellData[0]["name"].astype(str), "MT"))[0]
    mtGeneCounts = cellData[2][mtGeneNumbers, :].toarray().sum(axis=0)
    totals = np.asarray(cellData[2].sum(axis=0))[0]
    unique = np.asarray((cellData[2] > 0).sum(axis=0))[0]
    cellBarcodes =
    x = np.asarray(featureValues + [mtGeneCounts, totals, unique])
    d = pandas.DataFrame(data=x.transpose(),
                         columns=features + ["Mitochondrial UMIs", "Total UMIs", "Unique Features"])
    d["Percent MT"] = d["Mitochondrial UMIs"] / d["Total UMIs"]
    d["ScreenSeqTotal"] = d[features].sum(axis=1)
    d[[feature + "_norm" for feature in features]] = d[features].div(d["ScreenSeqTotal"], axis=0)
    passing = d[(d["Unique Features"] > 4000) & (d["Percent MT"] < 0.1) & (d["Total UMIs"] > 18000) & (
            d["Total UMIs"] < 80000)]


Analyze()
