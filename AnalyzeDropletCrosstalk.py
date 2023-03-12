import nd2
import numpy as np
import skimage.measure
import skimage.morphology
from PIL import Image
import pandas
from pathlib import Path

dataDirectory = Path(r"C:\Users\jonoj\Downloads\DropletCrosstalk\DropletCrosstalk")
paths = list(dataDirectory.iterdir())

areas = []
eccentricities = []
intensities = []
times = []

for time, path in enumerate(paths):
    nd2Data = nd2.ND2File(str(path.absolute())).asarray()
    for xy in range(nd2Data.shape[0]):
        print("t=%d, xy=%d" % (time, xy))
        bf = nd2Data[xy, 0]
        fitc = nd2Data[xy, 1]
        filteredBF = skimage.morphology.binary_opening(bf > 6000, skimage.morphology.disk(5))
        regions = skimage.measure.regionprops(skimage.measure.label(filteredBF))
        for region in regions:
            rows, cols = list(zip(*region.coords))
            if region.eccentricity > 0.75 or region.area < 2000 or region.area > 4000:
                filteredBF[rows, cols] = 0
            areas.append(region.area)
            eccentricities.append(region.eccentricity)
            intensities.append(np.sum(fitc[rows, cols] / 256))
            times.append(time)
        # bf = ((bf - np.min(bf)) / (np.max(bf) - np.min(bf))) * 255
        # bf = bf.astype(np.uint8)
        # overlay = np.stack([np.where(filteredBF, 255, bf),
        #                     np.where(filteredBF, 255, bf),
        #                     bf], axis=-1)
        # Image.fromarray(overlay).show()

data = pandas.DataFrame(np.asarray(list(zip(areas, intensities, eccentricities, times))),
                        columns=["Area", "Intensity", "Eccentricity", "Time"])
data.to_csv("CrossTalkData.csv", index=False)
