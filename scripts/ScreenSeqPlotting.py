import numpy as np
from scripts.ScreenSeq import ScreenSeqConfiguration
import math
import matplotlib.pyplot as plt
import seaborn
import scanpy as sc


def PlotQC(scData: sc.AnnData):
    axs = sc.pl.violin(scData, ["pct_counts_mt", "total_counts", "n_genes_by_counts"],
                       jitter=0.4, multi_panel=True, log=True, show=False)
    axs = axs.axes[0]
    axs[0].set_title("Percent mitochondrial")
    axs[0].set_xticks([])
    axs[0].set_ylabel("")
    axs[1].set_title("Transcript count")
    axs[1].set_xticks([])
    axs[1].set_ylabel("")
    axs[2].set_title("Unique genes")
    axs[2].set_xticks([])
    axs[2].set_ylabel("")

    plt.show()


def NormalizeKDEPlotHeight(ax):
    for collection in ax.collections:
        yCoords = collection.get_paths()[0].vertices[:, 1]
        yCoords = (yCoords - np.min(yCoords)) / (np.max(yCoords) - np.min(yCoords))
        collection.get_paths()[0].vertices[:, 1] = yCoords
    for line in ax.lines:
        yCoords = line.get_ydata()
        yCoords = (yCoords - np.min(yCoords)) / (np.max(yCoords) - np.min(yCoords))
        line.set_ydata(yCoords)


def PlotCalls(scData, config: ScreenSeqConfiguration):
    n = len(config.conditionBarcodeMapping.keys())
    rows = math.ceil(n / 4)
    fig, axs = plt.subplots(rows, 4, sharex="col", sharey="row",
                            gridspec_kw=dict(wspace=0, hspace=0))

    maxX = 0
    for (barcode, name), ax, color in zip(config.conditionBarcodeMapping.items(),
                                          axs.flatten(),
                                          seaborn.color_palette("husl", n_colors=n)):
        plt.sca(ax)
        barcodeFractions = scData.obs[barcode + "_frac"]
        barcodeCalls = scData.obs[barcode + "_call"]
        barcodeCalls = barcodeCalls[barcodeFractions > 0]
        barcodeFractions = barcodeFractions[barcodeFractions > 0]
        kdeAxis = seaborn.kdeplot(barcodeFractions, ax=ax, fill=False,
                                  color=[0, 0, 0, 0.2], alpha=1, clip=[0, 1])
        NormalizeKDEPlotHeight(kdeAxis)

        vertices = np.stack(kdeAxis.lines[0].get_data(), axis=-1)
        verticesClustered = scData.uns["clusterers"][barcode].predict(vertices[:, 0:1])
        splitCoords = [0] + list(np.where(np.diff(verticesClustered) != 0)[0] + 1) + [-1]
        segments = zip(splitCoords[0:], splitCoords[1:])
        for start, end in segments:
            end = min(-1, end + 1)
            pts = vertices[start:end, :]
            ax.fill_between(x=pts[:, 0], y1=pts[:, 1],
                            color=[0.7, 0.7, 0.7] if verticesClustered[start] == 0 else color,
                            alpha=1)

        plt.scatter(x=barcodeFractions[barcodeCalls == False],
                    y=np.random.uniform(-1, 1, np.count_nonzero(barcodeCalls == False)) * 0.025 - 0.05,
                    color=[0.7, 0.7, 0.7], s=0.25)
        plt.scatter(x=barcodeFractions[barcodeCalls == True],
                    y=np.random.uniform(-1, 1, np.count_nonzero(barcodeCalls == True)) * 0.025 - 0.05,
                    color=color, s=0.25)
        ax.set_yticks([])
        ax.set_ylabel(None)
        ax.set_xlabel(None)
        ax.set_ylim([-0.1, 1.1])
        maxX = max(maxX, barcodeFractions.max())
        ax.set_xlim([0, maxX*1.1])
        plt.text(np.mean(ax.get_xlim()), 0.9, name, ha="center", color=color)
    fig.supxlabel("Relative abundance")
    fig.supylabel("Normalized cell counts")
    plt.tight_layout()
    plt.show()