import fcsparser
import pandas
import seaborn
import matplotlib.pyplot as plt
import numpy as np


def GetData(path: str):
    return fcsparser.parse(path, reformat_meta=True)[1]


cd298 = GetData(r"Data\export_K562_CD298_Cells.fcs")
isotype = GetData(r"Data\export_K562_Isotype_Cells.fcs")
d = "640nm_670_30_APC-A"

cd298 = np.asarray(cd298[d])
isotype = np.asarray(isotype[d])

df = pandas.DataFrame(cd298, columns=["Cy5 intensity"])
df["Condition"] = "anti-CD298-Cy5"
df["Weight"] = 0.8

df2 = pandas.DataFrame(isotype, columns=["Cy5 intensity"])
df2["Condition"] = "Isotype-Cy5"
df2["Weight"] = 1
df = pandas.concat([df, df2])
df = df[df["Cy5 intensity"] > 1e-5]
seaborn.kdeplot(data=df, x="Cy5 intensity", hue="Condition", weights="Weight", fill=True, log_scale=True,
                palette=["green", "black"])

plt.yticks([])
plt.xlim([1, 1e6])
plt.title("Flow cytometry of K562 cell line")
plt.xlabel("Cy5 intensity")
plt.ylabel("Scaled density")
plt.show()
