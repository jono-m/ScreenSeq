import pickle
from typing import *
from scripts.DetectDropletsInVideo import Droplet
import matplotlib.pyplot as plt
import numpy as np

with open("droplets.pkl", "rb") as file:
    droplets: List[Droplet] = pickle.load(file)

colors = np.stack([droplet.meanColor for droplet in droplets])[None, :, :]
print(droplets[0].meanColor)
plt.imshow(colors / 255,
          aspect="auto")
plt.show()

