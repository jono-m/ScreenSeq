import numpy as np
import seaborn
import pandas
import matplotlib.pyplot as plt
import matplotlib.collections

data = pandas.read_csv("data.csv")

crosstalkData = data.loc[data["File"] == "Image20h001"]

crosstalkData = crosstalkData.loc[:, ("Droplet ID", "Channel", "Intensity")]
crosstalkData = crosstalkData[crosstalkData["Intensity"] < 25000]
crosstalkData = crosstalkData.pivot(index="Droplet ID", columns="Channel", values="Intensity").reset_index()

crosstalkData = crosstalkData.rename(columns={"4X_FITC_15%": "Fluorescein", "Rhodamine B": "Resorufin"})

crosstalkData["FluoresceinNorm"] = (crosstalkData["Fluorescein"] - np.min(crosstalkData["Fluorescein"])) / \
                                   (np.max(crosstalkData["Fluorescein"]) - np.min(
                                   crosstalkData["Fluorescein"]))
crosstalkData["ResorufinNorm"] = (crosstalkData["Resorufin"] - np.min(crosstalkData["Resorufin"])) / \
                                 (np.max(crosstalkData["Resorufin"]) - np.min(crosstalkData["Resorufin"]))

color_fluor = [0, 255, 26]
color_res = [255, 128, 0]

color = np.outer(crosstalkData["FluoresceinNorm"], color_fluor) + \
        np.outer(crosstalkData["ResorufinNorm"], color_res)
color = color / 255
color = color + 0.2
color = np.clip(color, 0, 1)

ax = seaborn.jointplot(data=crosstalkData, x="FluoresceinNorm", y="ResorufinNorm", kind="kde", zorder=0,
                       color=[0, 0, 0, 0.2])
plt.scatter(crosstalkData["FluoresceinNorm"], crosstalkData["ResorufinNorm"], c=color, marker="o",
            edgecolors=[0, 0, 0, 0.1], s=20, linewidths=0.5, zorder=1)

palette = seaborn.cubehelix_palette(start=1.8, hue=1, rot=0, dark=0.1, light=0.8, reverse=True,
                                    as_cmap=True)
points = ax.ax_marg_x.lines[0].get_path().vertices
ax.ax_marg_x.lines[0].remove()
segments = np.stack([points[:-1], points[1:]], axis=1)
lc = matplotlib.collections.LineCollection(segments,
                                           cmap=matplotlib.colors.LinearSegmentedColormap("Fluor", {
                                               "red": [(0, 0, 0),
                                                       (1, 0, 0)],
                                               "green": [(0, 0, 0),
                                                         (1, 1, 1)],
                                               "blue": [(0, 0, 0),
                                                        (1, 0.1, 0.1)]
                                           }))
lc.set_array(points[:, 0])
lc.set_linewidth(2)
ax.ax_marg_x.add_collection(lc)

points = ax.ax_marg_y.lines[0].get_path().vertices
ax.ax_marg_y.lines[0].remove()
segments = np.stack([points[:-1], points[1:]], axis=1)
lc = matplotlib.collections.LineCollection(segments,
                                           cmap=matplotlib.colors.LinearSegmentedColormap("Fluor", {
                                               "red": [(0, 0, 0),
                                                       (1, 1, 1)],
                                               "green": [(0, 0, 0),
                                                         (1, 0.5, 0.5)],
                                               "blue": [(0, 0, 0),
                                                        (1, 0, 0)]
                                           }))
lc.set_array(points[:, 1])
lc.set_linewidth(2)
ax.ax_marg_y.add_collection(lc)

plt.xlabel("Droplet Fluorescein Intensity (au)")
plt.ylabel("Droplet Resorufin Intensity (au)")
plt.show()
