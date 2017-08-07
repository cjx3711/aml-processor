from manifestExtraction import *
from unpairedFastqProc import *
from ReadPairer import *

fq1 = "MICROTEST_AD01_S1_L001_R1_001.fastq"
fq2 = "MICROTEST_AD01_S1_L001_R2_001.fastq"
inDir = "data/Raw/"
outDir = "data/Processed/"

numThreads = 12
readPairer = ReadPairer(probabilistic = False)


def pairToJ3X(fq1, fq2, inDir, outDir):
    with open(inDir + fq1) as fq1File, open(inDir + fq2) as fq2File:
        with open(outDir + fq1[:(fq1.find("_R1_"))] + "PAIRED.j3x", "w+", newline = "") as outFile:
            # Creates iterators which deliver the 4 lines of each FASTQ read as a zip (ID, Sequence, Blank, Quality)
            fq1Iter, fq2Iter = grouper(fq1File, 4), grouper(fq2File, 4)
            with ProcessPoolExecutor(numThreads) as processManager:
                # Calls alignAndMerge(FASTQ1's (ID, Sequence, Blank, Quality), FASTQ2's (ID, Sequence, Blank, Quality))
                for x in processManager.map(readPairer.alignAndMerge, fq1Iter, fq2Iter, chunksize = 250):
                    outFile.write(x)
                    outFile.write("\n\n")
                outFile.close()


if __name__ ==  "__main__": pairToJ3X(fq1, fq2, inDir, outDir)
