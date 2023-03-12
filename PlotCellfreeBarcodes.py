import pandas as pd
import matplotlib.pyplot as plt
import seaborn

data = pd.read_csv(
    r"C:\Users\jonoj\Repositories\ScreenSeq\Data\JMM_ScreenSeq_CellFree_10_7_2022.txt")
data = data.set_index(["Sample", "SamplePos", "Cycle"])
meanWater = data.loc["Water"].groupby("Cycle").mean()
data["SS1 Probe Intensity"] = data["465-510"] - meanWater["465-510"]
data["SS2 Probe Intensity"] = data["618-660"] - meanWater["618-660"]
data = data.drop("Water")

plt.subplot(1, 2, 1)
seaborn.lineplot(data, x="Cycle", y="SS1 Probe Intensity", hue="Sample", err_style="bars")

plt.subplot(1, 2, 2)
seaborn.lineplot(data, x="Cycle", y="SS2 Probe Intensity", hue="Sample", err_style="bars")

plt.show()