from manifestExtraction import *
from unpairedFastqProc import *
from pairedFastqProc import *

fileName = "AD01_S1_L001_R1_001.fastq"
fq1 = "MINITEST_AD01_S1_L001_R1_001.fastq"
fq2 = "MINITEST_AD01_S1_L001_R2_001.fastq"
inDir = "Reference Data/Raw/"
outDir = "Reference Data/Processed/"

genSignedReads(inDir, outDir)

#if __name__ ==  "__main__": pairToJ3X(fq1, fq2, inDir, outDir)
#if __name__ ==  "__main__": fastqToAny(fileName, inDir, outDir, "j3x")