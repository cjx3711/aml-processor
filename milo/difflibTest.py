from j3xUtils import *

def j3xIDAmp(j3x, inDir, outDir):
        with open(inDir + j3x) as j3xFile:
            with open(outDir + j3x + "IDed.j3x", "w+", newline = "") as outFile:
                # Creates iterator which deliver the 4 lines of each j3x read as a zip (ID, Sequence, Quality, Blank)
                j3xIter = grouper(j3xFile, 4)
                with ProcessPoolExecutor(numThreads) as processManager:
                    for x in processManager.map(idMe, j3xIter, chunksize = 250):
                        outFile.write(x)
                        outFile.write("\n\n")
                    outFile.close()