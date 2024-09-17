import numpy as np
import pandas
import matplotlib.pyplot as plt
import seaborn
import pathlib
import string

# Import Excel and format as dataframe.
data = pandas.read_excel(pathlib.Path(r"SS_Plate_3_18_2023.xlsx"), skiprows=51, usecols="C:M",
                         nrows=21,
                         header=None).to_numpy().flatten()
data = pandas.DataFrame(data, columns=["Luminescence"],
                        index=pandas.MultiIndex.from_product(
                            [list(string.ascii_uppercase[:7]), range(3), range(11)],
                            names=["Row", "Plate", "Column"])).reset_index()

# Conditions for each column:
conditions = pandas.DataFrame(data["Column"].map(
    lambda column: [["Tanespimycin", "vDMSO", 0.1],
                    ["Cytarabine", "vPBS", 1],
                    ["Dasatinib", "vDMSO", 0.1],
                    ["Homoharringtonine", "vDMSO", 0.1],
                    ["Hydroxyurea", "vPBS", 1],
                    ["Imatinib", "vDMSO", 1],
                    ["Ruxolitinib", "vDMSO", 0.1],
                    ["Vorinostat", "vDMSO", 1],
                    ["vDMSO", "N/A", np.nan],
                    ["vPBS", "N/A", np.nan],
                    ["+", "N/A", np.nan]][column]).tolist())
data[["Drug", "Vehicle", "Concentration"]] = conditions

# Serial dilution over rows. Positive control is not diluted.
data["Dilution"] = data["Row"].map(
    lambda row: np.power(0.1, string.ascii_uppercase.index(row)))
data["Dilution"] = np.where(data["Drug"] == "+", 1, data["Dilution"])
data["Concentration"] = data["Concentration"] * data["Dilution"] * 1000 * 0.1

positiveControl = data.query("Drug == '+'").groupby("Plate")["Luminescence"].mean()
negativeControl = data.query("Drug == 'vDMSO' | Drug == 'vPBS'").groupby(
    ["Plate", "Drug", "Dilution"])["Luminescence"].mean()

data = data.query("Vehicle != 'N/A'").reset_index()
data["Neg"] = list(negativeControl[zip(data["Plate"], data["Vehicle"], data["Dilution"])])
data["Pos"] = list(positiveControl[data["Plate"]])
data["Viability"] = (data["Luminescence"] - data["Pos"]) / (data["Neg"] - data["Pos"])

seaborn.lineplot(data=data, x="Concentration", y="Viability", hue="Drug")
plt.xscale("log")
plt.xlabel(r"Concentration ($\mu$g/mL)")
plt.show()
