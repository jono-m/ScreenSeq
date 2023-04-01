import numpy as np

from Utilities import *
from nd2reader import ND2Reader
import seaborn
import pandas
import matplotlib.pyplot as plt

dropletsAll = pandas.DataFrame()
for i in range(1, 9):
    images = []
    with ND2Reader(r"W:\Jono\FluorImages\Constant" + str(i) + ".nd2") as nd2file:
        for v in range(nd2file.sizes['v']):
            images.append((nd2file.get_frame_2D(v=v, c=0),
                           nd2file.get_frame_2D(v=v, c=1)))

    droplets = []
    for image in images:
        for droplet in FindDroplets(image[0], 200, 500, 0.5, False):
            intensity = np.mean(image[1][droplet.coords_np])
            droplet.intensity = intensity
            droplets.append(droplet)
    data = np.asarray([(droplet.intensity, i) for droplet in droplets])
    droplets = pandas.DataFrame(data=data,
                                columns=["Intensity", "Input"])
    dropletsAll = pandas.concat([dropletsAll, droplets], axis=0)
    print(i)

dropletsAll.to_csv("data.csv")
print(dropletsAll)
