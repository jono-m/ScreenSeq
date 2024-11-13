import sys
import pathlib
import os
from DropletVideoProcessor import DropletVideoProcessor

dvp = DropletVideoProcessor()
dvp.Read(r"../data/ColorMicroscopy/dropletVideo.mp4",
              29.73)

for p in pathlib.Path("output").iterdir():
    if p.is_file():
        os.remove(p)

for i in [0, 100, 200, 300, 750]:

    image = dvp.PreprocessImage(-0.5, (0.1, 0.9, .3, .7), 1.2, dvp.images[i])
    image.save("output/image_" + str(i) + ".png")