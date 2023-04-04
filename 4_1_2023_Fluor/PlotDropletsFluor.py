import numpy as np
import seaborn
import pandas
import matplotlib.pyplot as plt

data = pandas.read_csv("data.csv")

data["Intensity"] = data["Intensity"] / data["Intensity"].max()
ax = seaborn.boxplot(data=data, x="Intensity", orient="h", y="Inlet Number")
ax.legend_ = None
plt.ylabel("Inlet")
plt.xlabel("Droplet Fluorescence Intensity (au)")
plt.xlim([0, 1])
plt.yticks([])
plt.show()
