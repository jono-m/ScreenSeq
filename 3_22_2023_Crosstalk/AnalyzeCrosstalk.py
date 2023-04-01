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


images_neg = []
images_pos = []
images_mix = []
with nd2reader.ND2Reader(r"C:\Users\jonoj\Documents\Crosstalk20h001.nd2") as nd2file:
    for i in range(12):
        l = images_neg if i < 4 else (images_mix if i < 8 else images_pos)
        l.append((nd2file.get_frame_2D(v=i, c=0),
                  nd2file.get_frame_2D(v=i, c=1)))

intensities_neg = []
intensities_pos = []
intensities_mix = []


def Extract(inL, outL):
    for fitc, bf in inL:
        thresh = skimage.measure.label(bf > skimage.filters.threshold_mean(bf))
        regions = skimage.measure.regionprops(thresh)
        goodRegions = [r for r in regions if 12 < r.equivalent_diameter_area < 16]
        for r in goodRegions:
            c = np.asarray(r.coords)
            intensity = fitc[tuple(c.T)]
            outL.append(np.mean(intensity))


[Extract(inL, outL) for inL, outL in zip([images_neg, images_pos, images_mix],
                                         [intensities_neg, intensities_pos, intensities_mix])]
minI = np.min(intensities_neg + intensities_pos + intensities_mix)
maxI = np.max(intensities_neg + intensities_pos + intensities_mix)
intensities_neg = 100*(np.asarray(intensities_neg) - minI) / (maxI - minI)+1
intensities_pos = 100*(np.asarray(intensities_pos) - minI) / (maxI - minI)+1
intensities_mix = 100*(np.asarray(intensities_mix) - minI) / (maxI - minI)+1

fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True)

# bins = list(np.linspace(0, 10000, 100)) + list(np.linspace(10000, np.max(intensities_pos), 100))
plotLog = False
bins = 10 ** np.linspace(0, 2, 20) if plotLog else np.linspace(1, 100, 100)
plt.sca(axes[0])

seaborn.histplot(intensities_neg, color="black", bins=bins, kde=True)
seaborn.histplot(intensities_pos, color="#00FF00", bins=bins, kde=True)
plt.legend(["PBS droplets", "Fluorescein droplets"])
plt.sca(axes[1])
seaborn.histplot(intensities_mix, color="#008800", bins=bins, kde=True)
plt.legend(["Mixed droplets (24 hr)"])
plt.ylabel("Droplet count")
plt.xlabel("Fluorescence Intensity (au)")
plt.ylabel("Droplet count")
axes[1].invert_yaxis()
plt.subplots_adjust(hspace=0)
if plotLog:
    plt.xscale("log")
plt.show()
