import numpy as np
import seaborn
import pandas
import matplotlib.pyplot as plt

generator = np.random.Generator(np.random.PCG64())


def Sample(cell):
    receptorCount = int(cell[0])
    barcodeProbabilities = cell[1:]
    barcodes = range(len(cell[1:]))
    receptorStates = generator.choice(a=barcodes, size=receptorCount, p=barcodeProbabilities)
    barcodeCounts = np.bincount(receptorStates, minlength=len(barcodes))
    return barcodeCounts


def Lerp(a, b, t):
    return (b - a) * t + a


def CollectDropletsCommonComp(numberOfCells, valveStates, existing=None):
    # Describe the cells present
    receptorsPerCell = generator.normal(receptorsPerCellMean, receptorsPerCellSD, numberOfCells)
    receptorsPerCell = np.expand_dims(np.round(np.clip(receptorsPerCell, a_min=0, a_max=None)), -1)

    numberOfBarcodes = len(valveStates)

    barcodeConcentrations = np.random.normal(
        [Lerp(barcodeOffMean, barcodeOnMean, x) for x in valveStates],
        [Lerp(barcodeOffSD, barcodeOnSD, x) for x in valveStates],
        [numberOfCells, numberOfBarcodes])
    barcodeConcentrations = np.clip(barcodeConcentrations, a_min=0, a_max=None)

    numberOfOffBarcodes = sum([(1 - x) for x in valveStates])
    compConcentrations = np.random.normal(
        barcodeOnMean * numberOfOffBarcodes,
        barcodeOnSD * numberOfBarcodes,
        [numberOfCells, 1])
    compConcentrations = np.clip(compConcentrations, a_min=0, a_max=None)

    # Each cell receptor samples from the distribution of barcodes in the droplet to determine what
    # it is bound to.
    concentrations = np.concatenate([barcodeConcentrations, compConcentrations], axis=1)
    probabilities = concentrations / np.expand_dims(concentrations.sum(axis=1), axis=-1)
    cells = np.concatenate([receptorsPerCell, probabilities], axis=1)
    barcodeCounts = np.apply_along_axis(Sample, 1, cells)

    if existing is not None:
        return np.concatenate([existing, barcodeCounts], axis=0)
    return barcodeCounts


def CollectDropletsUniqueComp(numberOfCells, valveStates, existing=None):
    # Describe the cells present
    receptorsPerCell = generator.normal(receptorsPerCellMean, receptorsPerCellSD, numberOfCells)
    receptorsPerCell = np.expand_dims(np.round(np.clip(receptorsPerCell, a_min=0, a_max=None)), -1)

    numberOfBarcodes = len(valveStates)

    barcodeConcentrations = np.random.normal(
        [Lerp(barcodeOffMean, barcodeOnMean, x) for x in valveStates],
        [Lerp(barcodeOffSD, barcodeOnSD, x) for x in valveStates],
        [numberOfCells, numberOfBarcodes])
    barcodeConcentrations = np.clip(barcodeConcentrations, a_min=0, a_max=None)

    compConcentrations = np.random.normal(
        [Lerp(barcodeOffMean, barcodeOnMean, 1 - x) for x in valveStates],
        [Lerp(barcodeOffSD, barcodeOnSD, 1 - x) for x in valveStates],
        [numberOfCells, numberOfBarcodes])
    compConcentrations = np.clip(compConcentrations, a_min=0, a_max=None)

    # Each cell receptor samples from the distribution of barcodes in the droplet to determine what
    # it is bound to.
    concentrations = np.concatenate([barcodeConcentrations, compConcentrations], axis=1)
    probabilities = concentrations / np.expand_dims(concentrations.sum(axis=1), axis=-1)
    cells = np.concatenate([receptorsPerCell, probabilities], axis=1)
    barcodeCounts = np.apply_along_axis(Sample, 1, cells)

    if existing is not None:
        return np.concatenate([existing, barcodeCounts], axis=0)
    return barcodeCounts


def Normalize(x):
    return x / np.sum(x, axis=1, keepdims=True)


def NormalizeComp(x):
    b = int(x.shape[1] / 2)
    xs = x[:, :b]
    cs = x[:, b:]
    return xs / (xs + cs)


def Run(collectFunc, normFunc):
    droplets = None
    droplets = collectFunc(1000, [1] + [0] * 9, existing=droplets)
    droplets = collectFunc(1000, [0.75] + [0] * 9, existing=droplets)
    droplets = collectFunc(1000, [0.5] + [0] * 9, existing=droplets)
    droplets = collectFunc(1000, [0.25] + [0] * 9, existing=droplets)
    droplets = collectFunc(1000, [0] + [0] * 9, existing=droplets)
    droplets = normFunc(droplets)
    seaborn.stripplot(droplets[:, 0], orient="h")


def Common():
    Run(CollectDropletsCommonComp, Normalize)


def Unique():
    Run(CollectDropletsUniqueComp, NormalizeComp)


def Test():
    r = np.linspace(0, 1, 100)
    for a in r:
        droplets = NormalizeComp(CollectDropletsUniqueComp(1000, [a] + [0] * 9, existing=None))
        seaborn.kdeplot(droplets[:, 0])

    plt.legend(["%.2f" % a for a in r])
    plt.show()


receptorsPerCellMean = 800
receptorsPerCellSD = 0

# Describe the Screen-seq barcodes
oligosPerAntibodyMean = 5
oligosPerAntibodySD = 0.5

barcodeOnMean = 10000
barcodeOnSD = 0
barcodeOffMean = 0
barcodeOffSD = 0

Test()