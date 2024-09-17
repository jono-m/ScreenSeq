import pathlib
import sys
import numpy as np

import PIL.Image
from PIL import Image, ImageFilter
import matplotlib.pyplot as plt
import pandas

imagePaths = list(pathlib.Path(r"Cropped/").iterdir())

colors = np.zeros([len(imagePaths), 8])

for i,imagePath in enumerate(imagePaths):
    c, m, y, k = imagePath.stem.split("_")[1:-1]
    colors[i, :4] = [c, m, y, k]

    imageRaw = Image.open(str(imagePath))
    image = np.asarray(imageRaw)
    white = image.max(axis=2)
    r = image[:, :, 0]
    g = image[:, :, 1]
    b = image[:, :, 2]
    cyan = 255 - r
    magenta = 255 - g
    yellow = 255 - b
    key = 1-white
    cumulativePigment = [cyan.sum(), magenta.sum(), yellow.sum(), key.sum()]

    colors[i, 4:] = [r.sum(), g.sum(), b.sum(), 0]

print("here")