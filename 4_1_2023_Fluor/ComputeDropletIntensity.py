import numpy as np

from nd2reader import ND2Reader
import pandas
from PIL import Image
import skimage
from pathlib import Path

data = []
nd2Directory = Path(r"W:\Jono\Fluorescence2")
nd2Paths = [path for path in nd2Directory.iterdir() if path.suffix == ".nd2"]
segmentedPaths = Path(r"BF\Segmented").iterdir()
dropletID = 0
for nd2Path, segmentationPath in zip(nd2Paths, segmentedPaths):
    print("Analyzing %s and %s" % (nd2Path.name, segmentationPath.name))
    with ND2Reader(str(nd2Path.absolute())) as nd2file:
        with Image.open(segmentationPath) as segmentedImage:
            channels = nd2file.metadata['channels']
            nd2file.iter_axes = 'c'
            regions = skimage.measure.regionprops(skimage.measure.label(np.asarray(segmentedImage)))
            channelImages = [nd2file.get_frame(i) for i, _ in enumerate(channels)]
            for region in regions:
                for channel, channelImage in zip(channels, channelImages):
                    c = np.asarray(region.coords)
                    intensity = np.mean(channelImage[tuple(c.T)])
                    data.append((dropletID, segmentationPath.stem, channel, intensity))
                dropletID += 1

pandas.DataFrame(data=np.asarray(data), columns=["Droplet ID", "File", "Channel", "Intensity"]).to_csv("data.csv")
