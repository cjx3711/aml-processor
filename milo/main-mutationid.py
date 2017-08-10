from manifestExtraction import *
from concurrent.futures import *
from MutationFinder import *
import json
from pprint import pprint
import time

inDir = "data/2-paired/"
outDir = "data/3-mutations/"

numThreads = 12
chunksize = 250
mutationFinder = MutationFinder()

def run():
    minMutationCount = 2
    with open('config.json') as config_file: 
        config_data = json.load(config_file)
        if ( 'minMutationCount' in config_data ):
            minMutationCount = config_data['minMutationCount']
        if ( 'chunksize' in config_data ):
            chunksize = config_data['chunksize']
    
    print("Milo Mutation Identifier")
    print("Minimum Mutation Count: {0}".format(minMutationCount))
    print("Chunksize (Process Pool): {0}".format(chunksize))
    print()    
    
    with open('files.json') as file_list_file:    
        filenameArray = json.load(file_list_file)

        for filenames in filenameArray:
            pairedFile, mutationFile = readFilenames(filenames)
            
            if (pairedFile == '' or mutationFile == ''):
                print('Please set the keys "paired" and "mutation" in the config.json file')
                return
        
        for filenames in filenameArray:
            pairedFile, mutationFile = readFilenames(filenames)
            
            # Don't understand the __name__ thing, but it's required according to SO
            if __name__ ==  "__main__": mutationID(pairedFile, mutationFile, inDir, outDir, minMutationCount)
            
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
            print("{0} Crunching {1}".format(time.strftime('%X %d %b %Y') ,pairedFile))
            start = time.time()
            mutationFinder.reinit()
            # Creates iterators which deliver the 4 lines of each FASTQ read as a zip (ID, Sequence, Blank, Quality)
            inFileIter = grouper(inFile, 4)
            with ProcessPoolExecutor(numThreads) as processManager:
                # Calls alignAndMerge(FASTQ1's (ID, Sequence, Blank, Quality), FASTQ2's (ID, Sequence, Blank, Quality))
                for ampliconID, mutationHash in processManager.map(mutationFinder.identifyMutations, inFileIter, chunksize = chunksize):
                    mutationFinder.putMutationMap(ampliconID, mutationHash)
                
                print("{0} Dumping {1}".format(time.strftime('%X %d %b %Y'), mutationFile))
                mutationList = mutationFinder.extractHighestOccuringMutations(minMutationCount)
            
                for mutationOccurence in mutationList:    
                    outFile.write(str(mutationOccurence[1]))
                    outFile.write(', ')
                    outFile.write(mutationOccurence[0])
                    outFile.write('\n')
                
                outFile.close()
                print("{0} Dumped {1}".format(time.strftime('%X %d %b %Y'), mutationFile))
                print("Took {0}s".format(time.time() - start))
                
                
run()