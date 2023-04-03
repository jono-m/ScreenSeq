import numpy as np
import seaborn
import pandas
import matplotlib.pyplot as plt

data = pandas.read_csv("data.csv")

seaborn.boxplot(data=data, y="Intensity", x="Input")
plt.show()