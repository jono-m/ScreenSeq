import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn
import pathlib

data = pd.read_excel(pathlib.Path(r"SS_Drug_Screen_2_22_2023.xlsx"), skiprows=48,
                     usecols="B:E", nrows=72)
positiveControl = data[(data["Well"].str[0] == "H") & (data["Well"].str[1] != "1")][
    "Plate 1"].mean()
data = data[(data["Well"].str[0] != "H") | (data["Well"].str[1] == "1")]
drugsAndMgPerML = {"17-AAG": 1,
                   "Cytarabine": 10,
                   "Dasatinib": 5,
                   "Homoharringtonine": 1,
                   "Hydroxyurea": 10,
                   "Imatinib": 10,
                   "Ruxolitinib": 1,
                   "SAHA": 10,
                   "Vehicle": 1}

drugsAndIC5072 = {"17-AAG": 0.0002,
                  "Cytarabine": 0.0017,
                  "Dasatinib": 0.000002,
                  "Homoharringtonine": 0.0000136,
                  "Hydroxyurea": 0.0076,
                  "Imatinib": 0.00023,
                  "Ruxolitinib": 0.020,
                  "SAHA": 0.00055}

columns = np.asarray(sorted(drugsAndMgPerML.keys()))
data["Drug"] = columns[data["Well"].str[1].astype(int) - 1]
data["Dilution"] = 10 ** (data["Well"].str[0].map(str.upper).map(ord) - 65)
data["Concentration"] = data["Drug"].map(drugsAndMgPerML) / data["Dilution"] * 1000

data.loc[data["Well"] == "H1", "Drug"] = "DMSO20"
data.loc[data["Well"] == "H1", "Dilution"] = 1
data.loc[data["Well"] == "H1", "Concentration"] = 1

data = pd.wide_to_long(data, "Plate", sep=" ", i="Well", j="Replicate")
data = data.rename(columns={"Plate": "Luminescence"})

vehicle = data.groupby(["Drug", "Dilution"])["Luminescence"].mean()["Vehicle"]
# positiveControl = data.groupby(["Drug"])["Luminescence"].mean()["DMSO20"]

data["Viability"] = (data["Luminescence"] - positiveControl) / (
        data["Dilution"].map(vehicle) - positiveControl)

fig, axs = plt.subplots(2, 4)

for (drug, ax, color) in zip(columns, axs.flatten(), seaborn.color_palette()):
    plt.sca(ax)
    seaborn.lineplot(data=data.query("Drug == '%s'" % drug), x="Concentration", y="Viability",
                     color=color)
    seaborn.scatterplot(data=data.query("Drug == '%s'" % drug), x="Concentration", y="Viability",
                        color=color)
    plt.xlabel(r"Concentration ($\mu$g/mL)")
    plt.legend([drug])
    plt.ylim([-0.5, 1.5])
    plt.xscale("log")

fig.suptitle("K562 dose-response")
plt.show()
