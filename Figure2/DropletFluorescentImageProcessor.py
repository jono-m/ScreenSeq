import numpy as np
import skimage

def ComputeDropletIntensities(fluor, bf):
    droplets, segmented = _SegmentBF(bf)

    coords = [np.asarray(droplet.coords) for droplet in droplets]
    intensities = [fluor[coord[:, 0], coord[:, 1]].mean(axis=0) for coord in coords]
    return intensities


def _SegmentBF(grayscale):
    print("Segmenting...")
    grayscale = skimage.filters.gaussian(grayscale, 3)
    binary = grayscale > skimage.filters.threshold_local(grayscale, block_size=101)
    binary = skimage.segmentation.clear_border(binary)

    labeled = skimage.measure.label(binary)
    rps = skimage.measure.regionprops(labeled)

    def get_outliers(features):
        features = np.asarray(features)
        ix = np.where(abs(features - features.mean()) >= 1 * features.std())[0]
        return [rps[i].label for i in ix]

    outliers = get_outliers([rp.area for rp in rps])
    outliers += get_outliers([rp.solidity for rp in rps])

    labeled = np.where(np.isin(labeled, outliers), 0, labeled)

    segmented = labeled > 0

    return skimage.measure.regionprops(labeled), segmented
