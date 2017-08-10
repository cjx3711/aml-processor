from manifestExtraction import grouper
from MutationFinder import *
import json
from pprint import pprint
import time
import os
from multiprocessing import Pool, cpu_count
from tqdm import tqdm

inDir = "data/2-paired/"
outDir = "data/3-mutations/"

numThreads = cpu_count()
chunksize = 250
mutationFinder = MutationFinder()
bytesPerRead = 500 # Estimated

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
        filesize = os.path.getsize(inDir + pairedFile)
        estimatedReads = int(filesize / bytesPerRead)
        with open(outDir + mutationFile, "w+", newline = "") as outFile:
            print("{0} Crunching {1}".format(time.strftime('%X %d %b %Y') ,pairedFile))
            start = time.time()
            mutationFinder.reinit()
            # Creates iterators which deliver the 4 lines of each FASTQ read as a zip (ID, Sequence, Blank, Quality)
            inFileIter = grouper(inFile, 4)
            pool = Pool(numThreads)
            with tqdm(total=estimatedReads) as pbar:
                result = pool.imap_unordered(mutationFinder.identifyMutations, inFileIter, chunksize = chunksize)
                for i, data in tqdm(enumerate(result)):
                    ampliconID = data[0]
                    mutationHash = data[1]
                    referenceCoordinate = data[2]
                    mutationFinder.putMutationHash(ampliconID, mutationHash, referenceCoordinate)
                    pbar.update()

            pbar.close()
            pool.close()
            pool.join()
            
            print("{0} Dumping {1}".format(time.strftime('%X %d %b %Y'), mutationFile))
            mutationList = mutationFinder.extractHighestOccuringMutations(minMutationCount)
        
            for mutationOccurence in mutationList:
                mutationCount = mutationOccurence[1]
                mutationHash = mutationOccurence[0]
                ampliconID = int(mutationHash[:mutationHash.find(' ')])
                referenceAmplicon = mutationFinder.getReferenceAmplicon(ampliconID)
                mutationCoordinates = referenceAmplicon[1]
                appearanceCount = referenceAmplicon[2]
                outFile.write(str(mutationCount))
                outFile.write(' / ')
                outFile.write(str(appearanceCount))
                outFile.write('\t, ')
                outFile.write(mutationHash)
                outFile.write(', ')
                outFile.write(convertHashPositionsToCoordinates(mutationHash, mutationCoordinates))

                outFile.write('\n')
            
            outFile.close()
            print("{0} Dumped {1}".format(time.strftime('%X %d %b %Y'), mutationFile))
            print("Took {0}s".format(time.time() - start))
                
                
run()