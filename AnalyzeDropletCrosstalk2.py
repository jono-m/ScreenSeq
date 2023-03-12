import nd2
import numpy as np
import skimage.measure
import skimage.morphology
from PIL import Image
import pandas
from pathlib import Path
import math

dataPath = r"C:\Users\jonoj\Documents\Droplets12HR.nd2"
nd2Data = nd2.ND2File(dataPath).asarray()

sampleSet = ["PBS", "Fluorescein", "Mix"]

csvFile = open("CrossTalkData2.csv", "w+")
csvFile.write("Area,Intensity,Eccentricity,Sample\n")

for xy in range(nd2Data.shape[0]):
    sample = sampleSet[int(xy / 3)]
    bf = nd2Data[xy, 0]
    fitc = nd2Data[xy, 1]
    filteredBF = skimage.morphology.binary_opening(bf > 27000, skimage.morphology.disk(5))
    regions = skimage.measure.regionprops(skimage.measure.label(filteredBF))
    for region in regions:
        rows, cols = list(zip(*region.coords))
        if region.eccentricity > 0.8:
            filteredBF[rows, cols] = 0
        intensity = int(np.sum(fitc[rows, cols] / 256))
        csvFile.write("%d,%d,%f,%s\n" % (
            region.area, intensity, region.eccentricity, sample))
    #
    # bf = ((bf - np.min(bf)) / (np.max(bf) - np.min(bf))) * 255
    # bf = bf.astype(np.uint8)
    # overlay = np.stack([np.where(filteredBF, 255, bf),
    #                     np.where(filteredBF, 255, bf),
    #                     bf], axis=-1)
    # Image.fromarray(overlay).show()

csvFile.close()