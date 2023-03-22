import nd2reader
import numpy as np
import matplotlib.pyplot as plt

intensities = []
with nd2reader.ND2Reader(r"C:\Users\jonoj\Documents\dropletdata\ChangeCurve16Hz.nd2") as nd2file:
    for frame in nd2file.metadata['frames']:
        intensity = np.mean(nd2file.get_frame(frame))
        intensities.append(intensity)

print(len(intensities))

tPerFrame = 30 / len(intensities)
starts = [0, 70, 120, 180, 300]
freqs = [1, 2, 4, 8, 16]
cmap = plt.get_cmap('copper')
f, axs = plt.subplots(1, 5, sharey='all', sharex='all', gridspec_kw={'wspace': 0, 'hspace': 0})
for i, a in enumerate(axs):
    start = starts[i]
    end = start + 40
    plt.sca(a)
    xs = np.linspace(0, tPerFrame * 40, 40)
    plt.plot(xs, intensities[start:end], c=cmap(np.linspace(0, 1, 5)[i]))
    plt.yticks([])
    plt.xticks([0, 1, 2, 3])
    plt.title("%d Hz" % (freqs[i]))
plt.sca(axs[0])

f.add_subplot(111, frameon=False)
# hide tick and tick label of the big axes
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.grid(False)
plt.ylabel("Intensity (au)")
plt.xlabel("Time (sec)")
plt.show()
