import nd2reader
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
import pandas
import seaborn
import skimage


def normalizeImage(i):
    return (i - np.min(i)) / (np.max(i) - np.min(i))


def show(i):
    return Image.fromarray((normalizeImage(i) * 256).astype(int).astype(np.uint8)).show()


images = []
with nd2reader.ND2Reader(r"C:\Users\jonoj\Documents\MultiDroplets1.nd2") as nd2file:
    for i in range(4):
        images.append((nd2file.get_frame_2D(v=i, c=0),
                       nd2file.get_frame_2D(v=i, c=1)))

intensities = []


def Extract(inL, outL):
    for fitc, bf in inL:
        thresh = skimage.measure.label(bf > skimage.filters.threshold_mean(bf))
        regions = skimage.measure.regionprops(thresh)
        goodRegions = [r for r in regions if 10 < r.equivalent_diameter_area < 20]
        for r in goodRegions:
            c = np.asarray(r.coords)
            intensity = fitc[tuple(c.T)]
            outL.append(np.mean(intensity))


Extract(images, intensities)
# minI = np.min(intensities)
# maxI = np.max(intensities)
# intensities = 100 * (np.asarray(intensities) - minI) / (maxI - minI) + 1

plt.hist(intensities, bins=100)
plt.show()
