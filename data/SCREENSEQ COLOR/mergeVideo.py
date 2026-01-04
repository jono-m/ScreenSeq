from pathlib import Path
from scripts.DetectDropletsInVideo import DetectDropletsInVideo
import pandas
import matplotlib.pyplot as plt

directory = Path(r"./")
paths = list(directory.glob("*.avi"))

micronsPerPixel = 2.38
droplets = DetectDropletsInVideo(paths, (0.3, 0.35, 0, 0.35),
                                 debugMode=False, parallel=False, minDropletSize=1000, maxEccentricity=0.5)

df = pandas.DataFrame([droplet.rp["equivalent_diameter_area"] for droplet in droplets])

print("Droplet size: " + str(df.mean()*micronsPerPixel) + " µM (±" + str(df.std() * micronsPerPixel) + " µM)")
print("Diameter coefficient of variation: " + str(df.std()/df.mean()))

plt.hist(df, color="limegreen", bins="auto")

plt.show()