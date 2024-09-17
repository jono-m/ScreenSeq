import numpy as np
import seaborn
import pandas
import matplotlib.pyplot as plt

data = pandas.read_csv("data.csv")

pwmData = data.loc[data["File"].isin(["PWM100", "PWM25", "PWM75_2", "PWM50"]) & (data[
                                                                                     "Channel"] == "4X_FITC_15%")].copy()
pwmData.loc[pwmData["File"] == "PWM75_2", "File"] = "PWM75"
pwmData["Duty Cycle"] = pwmData["File"].str[3:].astype(int)
pwmData = pwmData.loc[:, ("Droplet ID", "Duty Cycle", "Intensity")]

gb = pwmData.groupby("Duty Cycle")["Intensity"]
mean = gb.transform('mean')
std = gb.transform('std')
zscore = np.abs(pwmData["Intensity"] - mean) / std
pwmData["ZScore"] = zscore

pwmData = pwmData[pwmData["ZScore"] <= 3].copy()
mean100 = gb.mean()[100]
pwmData["Intensity"] = pwmData["Intensity"] * 100 / mean100


def NormalizeKDEPlotHeight(ax):
    for collection in ax.collections:
        yCoords = collection.get_paths()[0].vertices[:, 1]
        yCoords = (yCoords - np.min(yCoords)) / (np.max(yCoords) - np.min(yCoords))
        collection.get_paths()[0].vertices[:, 1] = yCoords
    for line in ax.lines:
        yCoords = line.get_ydata()
        yCoords = (yCoords - np.min(yCoords)) / (np.max(yCoords) - np.min(yCoords))
        line.set_ydata(yCoords)


def kdePlotLabelPeaks(ax, y, labels, colors, *args, **kwargs):
    peaks = []
    for collection in ax.collections:
        yCoords = collection.get_paths()[0].vertices[:, 1]
        xCoords = collection.get_paths()[0].vertices[:, 0]
        peaks.append(xCoords[np.argmax(yCoords)])
    for line in ax.lines:
        yCoords = line.get_ydata()
        xCoords = line.get_xdata()
        peaks.append(xCoords[np.argmax(yCoords)])

    for peak, label, color in zip(reversed(peaks), labels, colors):
        plt.text(peak, y, label, *args, **kwargs, color=color, ha="center")

    return peaks


palette = seaborn.cubehelix_palette(n_colors=4, start=1.8, hue=1, rot=0,
                                    dark=0.1,
                                    light=.7, reverse=True)
ax = seaborn.kdeplot(data=pwmData, x="Intensity", hue="Duty Cycle", legend=False, fill=True,
                     palette=palette, alpha=0.8)
NormalizeKDEPlotHeight(ax)
peaks = kdePlotLabelPeaks(ax, 1.05, ["25%", "50%", "75%", "100%"], palette)

center = (max(peaks) + min(peaks)) / 2

plt.text(center, 1.15, "PWM duty cycle", ha="center")
plt.plot([min(peaks) - 10, max(peaks) + 10], [1.12, 1.12], color=[0, 0, 0, 0.5], linewidth=1)
hi = np.max(pwmData["Intensity"])
lo = 0
plt.xlim([lo, hi])
plt.ylim([0, 1.25])
plt.xticks([0, 25, 50, 75, 100])
plt.xlabel("Droplet fluorescence intensity (au)")
plt.ylabel("Normalized droplet count")
plt.yticks([])

plt.show()
