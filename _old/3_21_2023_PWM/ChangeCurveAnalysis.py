import numpy as np
from PIL import Image
import matplotlib.pyplot as plt

intensities = []
image = Image.open(r"C:\Users\jonoj\ChangeCurve\ChangeCurve.ome.tif")
for i in range(image.n_frames):
    image.seek(i)
    intensity = np.mean(np.asarray(image))
    intensities.append(intensity)

tPerFrame = 10 / len(intensities)
plt.plot(intensities)
plt.ylabel("Intensity (au)")
plt.xlabel("Time (sec)")
plt.show()
