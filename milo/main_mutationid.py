"""
Reads through the j3x and processes mutations and translocations
"""

from manifestExtraction import grouper
from MutationFinder import *
import json
from pprint import pprint
import time
import os
from multiprocessing import Pool, cpu_count
from tqdm import tqdm


class MainMutationID:
    def __init__(self):
        self.numThreads = cpu_count()
        self.bytesPerRead = 500 # Estimated

    def readConfig(self, configFile = 'config.json'):
        self.chunksize = 250
        self.minMutationCount = 3
        self.minTranslocationCount = 3
        with open(configFile) as config_file:
            config_data = json.load(config_file)
            if ( 'j4x_readThresholdForOutput' in config_data ):
                self.minMutationCount = config_data['j4x_readThresholdForOutput']
                self.minTranslocationCount = config_data['j4x_readThresholdForOutput']
            if ( 'j4x_readThresholdForMutationOutput' in config_data ):
                self.minMutationCount = config_data['j4x_readThresholdForMutationOutput']
            if ( 'j4x_readThresholdForTranslocationOutput' in config_data ):
                self.minTranslocationCount = config_data['j4x_readThresholdForTranslocationOutput']
            if ( 'chunksize' in config_data ):
                self.chunksize = config_data['chunksize']
    
    def test(self, inDir, outDir, configFile, referenceFile, filenameArray):
        self.readConfig(configFile)
        self.inDir = inDir
        self.outDir = outDir
        self.mutationFinder = MutationFinder(referenceFile)
        self.processFiles(filenameArray)
    
    def run(self):
        self.inDir = "data/2-paired/"
        self.outDir = "data/3-mutations/"
        self.mutationFinder = MutationFinder()
                
        self.readConfig()
        print("MILo Mutation Identifier")
        print("Minimum Mutation Count: {0}".format(self.minMutationCount))
        print("Chunksize (Process Pool): {0}".format(self.chunksize))
        print("Number of Threads: {0}".format(self.numThreads))
        print()
        
        with open('files.json') as file_list_file:  
            filenameArray = json.load(file_list_file)
            self.processFiles(filenameArray)
            
    def processFiles(self, filenameArray):
        for filenames in filenameArray:
            pairedFile, mutationFile, skip = self.readFilenames(filenames)
            if skip:
                continue
            if (pairedFile == '' or mutationFile == ''):
                print('Please set the keys "paired" and "mutation" in the files.json file')
                return
        
        for filenames in filenameArray:
            pairedFile, mutationFile, skip = self.readFilenames(filenames)
            if skip:
                continue
            # Don't understand the __name__ thing, but it's required according to SO
            self.mutationID(pairedFile, mutationFile)
        
    def readFilenames(self, filenames):
        pairedFile = mutationFile = ''
        skip = False
        
        if ( 'paired' in filenames ):
            pairedFile = filenames['paired']
        if ( 'mutation' in filenames ):
            mutationFile = filenames['mutation']
        if ( 'skip' in filenames ):
            skip = filenames['skip']
        return pairedFile, mutationFile, skip
        
    def mutationID(self, pairedFile, mutationFile):
        if not os.path.exists(self.outDir):
            os.makedirs(self.outDir)
        with open(self.inDir + pairedFile) as inFile:
            filesize = os.path.getsize(self.inDir + pairedFile)
            estimatedReads = int(filesize / self.bytesPerRead)
            with open(self.outDir + mutationFile, "w+", newline = "") as outFile:
                print("{0} Crunching {1}".format(time.strftime('%X %d %b %Y') ,pairedFile))
                start = time.time()
                self.mutationFinder.reinit()
                # Creates iterators which deliver the 4 lines of each FASTQ read as a zip (ID, Sequence, Blank, Quality)
                inFileIter = grouper(inFile, 4)
                pool = Pool(self.numThreads)
                with tqdm(total=estimatedReads) as pbar:
                    result = pool.imap_unordered(self.mutationFinder.identifyMutations, inFileIter, chunksize = self.chunksize)
                    for i, data in tqdm(enumerate(result)):
                        if data[0] == 'M':
                            ampliconID = data[1]
                            mutationHash = data[2]
                            referenceCoordinate = data[3]
                            readCount = data[4]
                            if ( ampliconID != None and mutationHash != None and referenceCoordinate != None and readCount != None ):
                                self.mutationFinder.putMutationHash(ampliconID, mutationHash, referenceCoordinate, readCount)
                        else: # Translocation
                            ampID1 = data[1]
                            ampID2 = data[2]
                            matchingBlocks = data[3]
                            readCount = data[4]
                            if ( ampID1 != None and ampID2 != None and matchingBlocks != None and readCount != None ):
                                self.mutationFinder.putTranslocationHash(ampID1, ampID2, matchingBlocks, readCount)
                        pbar.update()

                pbar.close()
                pool.close()
                pool.join()
                
                print("{0} Dumping {1}".format(time.strftime('%X %d %b %Y'), mutationFile))
                mutationList = self.mutationFinder.extractHighestOccuringMutations(self.minMutationCount)
                translocationList = self.mutationFinder.extractHighestOccuringTranslocations(self.minTranslocationCount)
                
                outFile.write('Mutations\n')
                for mutationOccurence in mutationList:
                    mutationCount = mutationOccurence[1]
                    mutationHash = mutationOccurence[0]
                    ampliconID = int(float(mutationHash[:mutationHash.find(' ')]))
                    referenceAmplicon = self.mutationFinder.getReferenceAmplicon(ampliconID)
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
                    
                outFile.write('Translocations\n')
                
                for translocationOccurence in translocationList:
                    translocationCount = translocationOccurence[1]
                    translocationHash = translocationOccurence[0]
                    outFile.write(str(translocationCount))
                    outFile.write('\t, ')
                    outFile.write(translocationHash)
                    outFile.write('\n')
                
                    
                outFile.write('Reference Amplicon Stats\n')
                referenceAmplicons = self.mutationFinder.getReferenceAmpliconArray()
                outFile.write('Count: ' + str(self.mutationFinder.referenceCount) + '\n')
                outFile.write('ampID, Reads, Mutations, Translocations, Sequence\n')
                
                total = 0
                for i, refAmp in enumerate(referenceAmplicons):
                    total += refAmp[2]
                    outFile.write(str(i+1) + ', ')
                    outFile.write(str(refAmp[2]) + ', ')
                    outFile.write(str(refAmp[3]) + ', ')
                    outFile.write(str(refAmp[4]) + ', ')
                    outFile.write(refAmp[0] + '\n')
                outFile.write('Total: ' + str(total))
                
                
                outFile.close()
                print("{0} Dumped {1}".format(time.strftime('%X %d %b %Y'), mutationFile))
                print("Took {0}s\n\n".format(time.time() - start))
                    
if __name__ ==  "__main__":
    main = MainMutationID()
    main.run()