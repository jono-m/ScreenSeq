import numpy as np
import pandas
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

data = pandas.read_csv("Calls.csv", index_col=0)
data = data[[column for column in data.columns if column not in ["Omacetaxine", "Vorinostat"]]]


def ComputeCounts(data):
    nBarcodes = len(data.columns)
    nConditions = 2 ** nBarcodes
    conditionCounts = []
    for condition in range(nConditions):
        states = [bool(int(x)) for x in format(condition, '0%db' % nBarcodes)]
        conditionCounts.append(np.count_nonzero(np.all(data == states, axis=1)))
    return conditionCounts


# Do the plotting
fig, axes = plt.subplots()
ax = axes


def DrawConditionsMatrix(ax, drugs, horizontal):
    ax.set_aspect(1)
    cellSize = 1
    cellSpacing = 1
    nDrugs = len(drugs)
    nConditions = 2 ** nDrugs
    d1 = nDrugs * (cellSize + cellSpacing) + cellSpacing
    d2 = nConditions * (cellSize + cellSpacing) + cellSpacing
    ax.set_xlim(0, d1 if horizontal else d2)
    ax.set_ylim(0, d2 if horizontal else d1)

    ax.yaxis.set_visible(horizontal)
    ax.yaxis.set_visible(not horizontal)
    ax.xaxis.set_tick_params(length=0, labeltop=True, labelbottom=False)
    ax.yaxis.set_tick_params(length=0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    for conditionNumber in range(nConditions):
        states = [bool(int(x)) for x in format(conditionNumber, '0%db' % nDrugs)]
        for col, state in enumerate(states):
            xy = (col * (cellSize + cellSpacing) + cellSpacing,
                  conditionNumber * (cellSize + cellSpacing) + cellSpacing)
            ax.add_patch(
                Rectangle(xy=xy, width=cellSize, height=cellSize,
                          facecolor=[0.8, 0, 0] if state else "white",
                          edgecolor=[0.4, 0, 0] if state else [0.6, 0.6, 0.6]))
        count = conditionCounts[conditionNumber]
        (x, y) = (7 * (cellSize + cellSpacing) + cellSpacing,
                  conditionNumber * (cellSize + cellSpacing) + cellSpacing)
        ax.text(x, y, str(count), ha="left", va="center")

    ax.set_xticks(
        [col * (cellSize + cellSpacing) + cellSpacing + (cellSize / 2) for col in
         range(nDrugs + nDataCols)],
        labels=list(data.columns) + ["n"], rotation=90)


plt.show()
