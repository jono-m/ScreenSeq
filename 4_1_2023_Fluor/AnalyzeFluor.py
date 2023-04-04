import numpy as np

from PIL import Image
import skimage
from nd2reader import ND2Reader
import seaborn
import pandas
import matplotlib.pyplot as plt


def Scale01(i):
    i = i.astype(float)
    return (i - np.min(i)) / (np.max(i) - np.min(i))


def Show(i):
    return Image.fromarray(np.floor(Scale01(i) * 255).astype(np.uint8)).show()


def Analyze(ix):
    with ND2Reader(r"C:\Users\jonoj\Documents\Fluorescence2\Inlet%d.nd2" % ix) as nd2file:
        image = nd2file.get_frame(0)
        thresholded = image > skimage.filters.threshold_local(image, block_size=21)
        thresholded = skimage.morphology.binary_opening(thresholded, skimage.morphology.disk(5))

        regions = skimage.measure.regionprops(skimage.measure.label(thresholded))

        for region in regions:
            circularity = region.area * 4 * np.pi / (region.perimeter_crofton ** 2)
            region.circularity = circularity
            region.coords_np = tuple(np.asarray(region.coords).T)
            region.intensity = np.mean(image[region.coords_np])

        overlay = (Scale01(image) * 255).astype(np.uint8)
        overlay = np.stack([overlay] * 3, axis=-1)
        circularityCutoff = 0.1
        for region in regions:
            if region.circularity > circularityCutoff and 10 <= region.axis_minor_length <= 25:
                overlay[region.coords_np[0], region.coords_np[1], 1] = 255

        intensities = [region.intensity for region in regions if region.circularity > circularityCutoff]
        data = np.asarray([[ix] * len(intensities), intensities]).T
        Image.fromarray(overlay).show()
    return pandas.DataFrame(columns=["Inlet Number", "Intensity"], data=data)


Analyze(4)
