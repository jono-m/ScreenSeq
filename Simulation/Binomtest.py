import numpy as np
import seaborn
import matplotlib.pyplot as plt


def Sample(n, p):
    a = np.random.choice([0, 1], n, p=[p, 1 - p])
    c = np.bincount(a, minlength=2)
    return c


samplesToTake = 1000
for i, n in enumerate([100, 1000]):
    p = 0.1
    allSuccesses = []
    for sampleNumber in range(samplesToTake):
        successes, _ = Sample(n, p)
        allSuccesses.append(successes)
    plt.subplot(2, 1, i + 1)
    seaborn.histplot(allSuccesses)
    plt.xlim([0, p * n * 2])

plt.show()
