
import sys

from scripts import ScreenSeq, ScreenSeqPlotting

path = r"../data/Sequencing/screenSeq_10_2025.h5"
scData = ScreenSeq.Load10X(path)
scData = ScreenSeq.FilterCells(scData, minTranscripts=1000, minUniqueGenes=1000, maxMitochondrialPercent=10)

config = ScreenSeq.ScreenSeqConfiguration()
config.AddConditionBarcode("ScreenSeq1", "PMA")
config.AddConditionBarcode("ScreenSeq2", "Ionomycin")
config.AddConditionBarcode("ScreenSeq3", "Anti-CD3")
config.AddConditionBarcode("ScreenSeq4", "Anti-CD28")
config.AddConditionBarcode("ScreenSeq5", "IL2")
config.AddConditionBarcode("ScreenSeq6", "IL4")
config.AddConditionBarcode("ScreenSeq7", "IL15")
config.AddConditionBarcode("ScreenSeq8", "IL21")
config.AddNormalizationBarcode("ScreenSeqN")

scData = ScreenSeq.CallConditions(scData, config)
ScreenSeqPlotting.PlotCalls(scData, config)
scData = ScreenSeq.Normalize(scData)
scData = ScreenSeq.Cluster(scData)

sys.exit()