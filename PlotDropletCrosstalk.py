import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

data = pd.read_csv("CrossTalkData.csv")
# data["Time"] = data["Time"].astype(int).astype(str) + " hr"
data["Droplet Fluorescence"] = data["Intensity"] / data["Area"]
good_droplets = data.query("Eccentricity < 0.75 & Area > 2000 & Area < 4000")

sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
# Initialize the FacetGrid object
pal = sns.cubehelix_palette(10, start=2, rot=0, light=.7)
g = sns.FacetGrid(good_droplets, row="Time", hue="Time", aspect=10, height=1, palette=pal)

# Draw the densities in a few steps
g.map(sns.kdeplot, "Droplet Fluorescence",
      bw_adjust=.5, clip_on=False,
      fill=True, alpha=1, linewidth=1.5)
g.map(sns.kdeplot, "Droplet Fluorescence", clip_on=False, color="w", lw=2, bw_adjust=.5)

# passing color=None to refline() uses the hue mapping
g.refline(y=0, linewidth=2, linestyle="-", color=None, clip_on=False)


def label(x, color, label):
    ax = plt.gca()
    ax.text(0, .2, label + " hr", fontweight="bold", color=color,
            ha="left", va="center", transform=ax.transAxes)


g.map(label, "Droplet Fluorescence")
# Set the subplots to overlap
g.figure.subplots_adjust(hspace=-.25)

# Remove axes details that don't play well with overlap
g.set_titles("")
g.set(yticks=[], ylabel="")
g.despine(bottom=True, left=True)
plt.show()
