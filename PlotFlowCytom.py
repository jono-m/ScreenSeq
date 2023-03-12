import pandas as pd
import matplotlib.pyplot as plt
import seaborn

col = "640nm_670_30_APC-A"
isoType = pd.read_csv(
    r"C:\Users\jonoj\Google Drive\Research\Screen-seq\Data\export_Sample_IsotypeControl.csv")
CD298 = pd.read_csv(
    r"C:\Users\jonoj\Google Drive\Research\Screen-seq\Data\export_Sample_CD298.csv")
isoType["Antibody"] = "Isotype"
CD298["Antibody"] = "CD298"
data = pd.concat([isoType, CD298])
seaborn.histplot(data, x=col, hue="Antibody")
plt.xlabel("Fluorescence Intensity on K562 Cells")
plt.show()
