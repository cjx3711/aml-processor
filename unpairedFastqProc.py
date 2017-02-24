"""
All FastQ processing routines that do not involve paired reads.
"""

from j3xUtils import *
from concurrent.futures import *

formatsToFuncs = {"j3x" : fastqToJ3X}
numThreads = 8

def fastqToAny(fastqName, inDir, outDir, convertTo):
        with open(inDir + fastqName) as fastqFile:
            with open(outDir + fastqName[:-6] + "." + convertTo, "w+", newline = "") as outFile:
                with ProcessPoolExecutor(numThreads) as processManager:
                    for x in processManager.map(formatsToFuncs[convertTo], grouper(fastqFile, 4), chunksize = 1200):
                        outFile.write(x)
                    outFile.close()