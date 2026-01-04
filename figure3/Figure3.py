import sys

import scripts.ScreenSeq as ss

path = r"../data/Sequencing/screenSeq_3_2023.h5"
data = ss.Load10X(path)
data = ss.FilterCells(data, maxMitochondrialPercent=7, minTranscripts=1000, minUniqueGenes=1000)

config = ss.ScreenSeqConfiguration()
config.AddConditionBarcode("ScreenSeq4", "MonoA")
config.AddConditionBarcode("ScreenSeq5", "MonoB")
config.AddConditionBarcode("ScreenSeq6", "Pair-1")
config.AddConditionBarcode("ScreenSeq7", "Pair-2")
config.AddConditionBarcode("ScreenSeq8", "Trio-1")
config.AddConditionBarcode("ScreenSeq9", "Trio-2")
config.AddConditionBarcode("ScreenSeq10", "Trio-3")
config.AddNormalizationBarcode("ScreenSeq11")
config.AddNormalizationBarcode("ScreenSeqN")

data = ss.CallConditions(data, config)
sys.exit()