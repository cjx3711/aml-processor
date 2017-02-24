from extractManifest import *
from unpairedFastqProc import *

fileName = "AD01_S1_L001_R1_001MINITEST.fastq"
inDir = "Reference Data/Raw/"
outDir = "Reference Data/Processed/"

if __name__ ==  "__main__": fastqToAny(fileName, inDir, outDir, "j3x")