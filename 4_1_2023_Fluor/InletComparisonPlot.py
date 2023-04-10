import numpy as np
import seaborn
import pandas
import matplotlib.pyplot as plt

data = pandas.read_csv("data.csv")
data = data.loc[data["File"].isin(["Inlet%d" % i for i in range(1, 9)])]
data = data.loc[data["Channel"] == "4X_FITC_15%"]
data = data.loc[data["Intensity"] > 5000].copy()

grpData = data.groupby("File")["Intensity"]
mean = grpData.transform("mean")
std = grpData.transform("std")
data["ZScore"] = np.abs(data["Intensity"] - mean) / std

data = data[data["ZScore"] <= 3].copy()
means = data["Intensity"].mean()
data.loc[:, "Intensity"] = (data["Intensity"] - means)/means
seaborn.violinplot(data=data, y="File", x="Intensity", orient="h", color="Green", showfliers=False)
plt.axvline(color="black")
plt.xlabel("Droplet fluorescence deviation from mean (%)")
plt.ylabel("Inlet")
plt.xlim([-10, 10])
plt.yticks(plt.gca().get_yticks(), labels=[str(x) for x in range(1, 9)])
plt.show()
