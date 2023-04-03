import numpy as np

from Utilities import *
from nd2reader import ND2Reader
import seaborn
import pandas
import matplotlib.pyplot as plt

with ND2Reader(r"C:\Users\jonoj\Documents\Fluorescence2\Inlet1.nd2") as nd2file:
    print("here")