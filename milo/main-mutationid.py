from manifestExtraction import *
from concurrent.futures import *
from MutationFinder import *
import json
from pprint import pprint
import time

inDir = "data/Processed/"
outDir = "data/Processed/"

numThreads = 12
mutationFinder = MutationFinder()

def run():
    minMutationCount = 2
    with open('config.json') as config_file: 
        config_data = json.load(config_file)
        if ( 'minMutationCount' in config_data ):
            minMutationCount = config_data['minMutationCount']
            
    with open('files.json') as file_list_file:    
        filenameArray = json.load(file_list_file)

        for filenames in filenameArray:
            pairedFile, mutationFile = readFilenames(filenames)
            
            if (pairedFile == '' or mutationFile == ''):
                print('Please set the keys "paired" and "mutation" in the config.json file')
                return
        
        for filenames in filenameArray:
            pairedFile, mutationFile = readFilenames(filenames)
            start = time.time()
            # Don't understand the __name__ thing, but it's required according to SO
            if __name__ ==  "__main__": mutationID(pairedFile, mutationFile, inDir, outDir, minMutationCount)
            end = time.time()
            print("Took {0}s".format(end - start))

def readFilenames(filenames):
    pairedFile = mutationFile = ''
    if ( 'paired' in filenames ):
        pairedFile = filenames['paired']
    if ( 'mutation' in filenames ):
        mutationFile = filenames['mutation']
    return pairedFile, mutationFile
    
def mutationID(pairedFile, mutationFile, inDir, outDir, minMutationCount):
    with open(inDir + pairedFile) as inFile:
        with open(outDir + mutationFile, "w+", newline = "") as outFile:
            print("Crunching {0}".format(pairedFile))
            mutationFinder.reinit()
            # Creates iterators which deliver the 4 lines of each FASTQ read as a zip (ID, Sequence, Blank, Quality)
            inFileIter = grouper(inFile, 4)
            with ProcessPoolExecutor(numThreads) as processManager:
                # Calls alignAndMerge(FASTQ1's (ID, Sequence, Blank, Quality), FASTQ2's (ID, Sequence, Blank, Quality))
                for ampliconID, mutationHash in processManager.map(mutationFinder.identifyMutations, inFileIter, chunksize = 250):
                    mutationFinder.putMutationMap(ampliconID, mutationHash)
                
                print("Dumping {0}".format(mutationFile))
                mutationList = mutationFinder.extractHighestOccuringMutations(minMutationCount)
            
                for mutationOccurence in mutationList:    
                    outFile.write(str(mutationOccurence[1]))
                    outFile.write(', ')
                    outFile.write(mutationOccurence[0])
                    outFile.write('\n')
                
                outFile.close()
                print("Dumped {0}".format(mutationFile))
                
run()