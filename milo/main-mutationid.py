from manifestExtraction import *
from concurrent.futures import *

from MutationFinder import *


j3x = "MICROTEST_AD01_S1_L001PAIRED.j3x"
inDir = "data/Processed/"
outDir = "data/Processed/"

numThreads = 4


mutationFinder = MutationFinder()
def mutationID(j3x, inDir, outDir):
    with open(inDir + j3x) as j3xFile:
        with open(outDir + j3x + '_MUTATIONS.j4x', "w+", newline = "") as outFile:
            # Creates iterators which deliver the 4 lines of each FASTQ read as a zip (ID, Sequence, Blank, Quality)
            j3xIter = grouper(j3xFile, 4)
            with ProcessPoolExecutor(numThreads) as processManager:
                # Calls alignAndMerge(FASTQ1's (ID, Sequence, Blank, Quality), FASTQ2's (ID, Sequence, Blank, Quality))
                for ampliconID, mutationHash in processManager.map(mutationFinder.identifyMutations, j3xIter, chunksize = 250):
                    mutationFinder.putMutationMap(ampliconID, mutationHash)
                    outFile.write(str(ampliconID))
                    outFile.write("\n\n")
                
                outFile.close()
                print("Done")
                tuples = list(mutationFinder.getMutationMap().items())
                tuples.sort(key=lambda tup: tup[1]) 
                print(tuples)
                
if __name__ ==  "__main__": mutationID(j3x, inDir, outDir)
