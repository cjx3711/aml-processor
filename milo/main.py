from manifestExtraction import *
from unpairedFastqProc import *
from pairedFastqProc import *

fileName = "AD01_S1_L001_R1_001.fastq"
fq1 = "AD01_S1_L001_R1_001MICROTEST.fastq"
fq2 = "AD01_S1_L001_R2_001MICROTEST.fastq"
inDir = "Reference Data/Raw/"
outDir = "Reference Data/Processed/"

if __name__ ==  "__main__": pairToJ3X(fq1, fq2, inDir, outDir)
#if __name__ ==  "__main__": fastqToAny(fileName, inDir, outDir, "j3x")