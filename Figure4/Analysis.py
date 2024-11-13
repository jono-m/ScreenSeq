import sys

from src import ScreenSeq, ScreenSeqPlotting

path = r"../data/Sequencing/screenSeq_9_2024.h5.h5"
scData = ScreenSeq.Load10X(path)
scData = ScreenSeq.FilterCells(scData, minTranscripts=1000, minUniqueGenes=1000)

config = ScreenSeq.ScreenSeqConfiguration()
config.AddConditionBarcode("ScreenSeq1", "LPS")
config.AddConditionBarcode("ScreenSeq2", "IFN-gamma")
config.AddConditionBarcode("ScreenSeq3", "TGF-beta")
config.AddConditionBarcode("ScreenSeq4", "TNF-alpha")
config.AddConditionBarcode("ScreenSeq5", "PMA")
config.AddConditionBarcode("ScreenSeq6", "IL6")
config.AddConditionBarcode("ScreenSeq7", "IL1-beta")
config.AddConditionBarcode("ScreenSeq8", "Mock")
config.AddNormalizationBarcode("ScreenSeqN")

scData = ScreenSeq.CallConditions(scData, config)
ScreenSeqPlotting.PlotCalls(scData, config)
scData = ScreenSeq.Normalize(scData)
scData = ScreenSeq.Cluster(scData)

sys.exit()