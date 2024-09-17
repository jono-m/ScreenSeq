import importlib
import Figure2_plot
import Figure2_analysis


def Run(a: Figure2_analysis.Analysis):
    importlib.reload(Figure2_plot)
    Figure2_plot.Plot(a)