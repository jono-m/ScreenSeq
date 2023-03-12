import pandas as pd
import seaborn
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

data = pd.read_csv("CrossTalkData2.csv")
data["Droplet Fluorescence"] = data["Intensity"] / data["Area"]
good_droplets = data.query("Eccentricity < 0.5 & Area > 200 & Area < 400")
seaborn.kdeplot(good_droplets, x="Droplet Fluorescence", hue="Droplet Sample")
plt.title("12-Hour Droplet Crosstalk")
plt.show()
