import pandas
import seaborn
import matplotlib.pyplot as plt
import numpy as np
from sklearn.cluster import MeanShift
import tables


def NormalizeKDEPlotHeight(ax):
    for collection in ax.collections:
        yCoords = collection.get_paths()[0].vertices[:, 1]
        yCoords = (yCoords - np.min(yCoords)) / (np.max(yCoords) - np.min(yCoords))
        collection.get_paths()[0].vertices[:, 1] = yCoords
    for line in ax.lines:
        yCoords = line.get_ydata()
        yCoords = (yCoords - np.min(yCoords)) / (np.max(yCoords) - np.min(yCoords))
        line.set_ydata(yCoords)


def CallSSData(data):
    calls = []
    clusterers = []
    for drug in data.columns:
        print("Clustering " + drug)
        barcodeAmounts = data[drug]

        clusterer = MeanShift(bin_seeding=True, min_bin_freq=round(len(barcodeAmounts) / 10)).fit(
            np.asarray(barcodeAmounts).reshape(-1, 1))
        calls.append((drug, clusterer.labels_ > 0))
        clusterers.append(clusterer)

    df = pandas.DataFrame(data=np.asarray([response for drug, response in calls]).T,
                          index=data.index,
                          columns=[drug for drug, response in calls])
    return df, clusterers


def PlotCalls(data, calls, clusterers):
    fig, axs = plt.subplots(2, 4, sharex="col", sharey="row", gridspec_kw=dict(wspace=0, hspace=0))

    for drug, ax, color, clusterer in zip(data.columns, axs.flatten(),
                                          seaborn.color_palette("husl", n_colors=len(data.columns)),
                                          clusterers):
        print("Plotting " + drug)
        plt.sca(ax)

        barcodeAmounts = data[drug]
        kdeAxis = seaborn.kdeplot(barcodeAmounts, ax=ax, fill=False, color=[0, 0, 0, 0.2],
                                  alpha=0.2,
                                  clip=[0, 100])
        NormalizeKDEPlotHeight(kdeAxis)

        vertices = np.stack(kdeAxis.lines[0].get_data(), axis=-1)
        verticesClustered = clusterer.predict(vertices[:, 0:1])
        splitCoords = [0] + list(np.where(np.diff(verticesClustered) != 0)[0] + 1) + [-1]
        segments = zip(splitCoords[0:], splitCoords[1:])
        for start, end in segments:
            cluster = verticesClustered[start]
            end = min(-1, end + 1)
            pts = vertices[start:end, :]
            ax.fill_between(x=pts[:, 0], y1=pts[:, 1],
                            color=[0.7, 0.7, 0.7] if cluster == 0 else color,
                            alpha=0.5)

        plt.scatter(x=data.loc[~calls[drug], drug],
                    y=np.random.uniform(-1, 1, np.count_nonzero(~calls[drug])) * 0.025 - 0.05,
                    color=[0.7, 0.7, 0.7], s=0.25)
        plt.scatter(x=data.loc[calls[drug], drug],
                    y=np.random.uniform(-1, 1, np.count_nonzero(calls[drug])) * 0.025 - 0.05,
                    color=color, s=0.25)
        ax.set_yticks([])
        ax.set_ylabel(None)
        ax.set_xlabel(None)
        ax.set_ylim([-0.1, 1.1])
        ax.set_xlim([0, 40])
        ax.set_xticks([0, 10, 20, 30])
        plt.text(np.mean(ax.get_xlim()), 0.9, drug, ha="center", color=color)

    # Turn off axis lines and ticks of the big subplot
    fig.supxlabel("Drug barcode relative abundance (%)")
    fig.supylabel("Normalized cell counts")
    plt.tight_layout()
    plt.show()
