import ScreenSeq

path = r"PlateScreen_7_2024.h5"
scData = ScreenSeq.Load10X(path)

scData = ScreenSeq.FilterCells(scData, maxMitochondrialPercent=5.75,
                               minTranscripts=2000, maxTranscripts=55000,
                               minUniqueGenes=1000)

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

ScreenSeq.Normalize(scData)
ScreenSeq.CellCycleScore(scData)
ScreenSeq.RegressCellCycle(scData)

scData = ScreenSeq.Cluster(scData)

print("test")