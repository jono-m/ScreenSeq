import sys

from scripts import ScreenSeq, ScreenSeqPlotting

path = r"_old/2_10_2024_Screen/ssdata_2_10_2024.h5"
scData = ScreenSeq.Load10X(path)
scData = ScreenSeq.FilterCells(scData, minTranscripts=1000, minUniqueGenes=1000)
config = ScreenSeq.ScreenSeqConfiguration()
config.AddConditionBarcode("ScreenSeq4", "E. Coli")
config.AddConditionBarcode("ScreenSeq7", "IFN-gamma")
config.AddConditionBarcode("ScreenSeq8", "IFN-beta")
config.AddConditionBarcode("ScreenSeq9", "TNF-alpha")
config.AddConditionBarcode("ScreenSeq10", "IL1-beta")
config.AddNormalizationBarcode("ScreenSeq3")

scData = ScreenSeq.CallConditions(scData, config)
ScreenSeqPlotting.PlotCalls(scData, config)
scData = ScreenSeq.Normalize(scData)
scData = ScreenSeq.Cluster(scData)

sys.exit()