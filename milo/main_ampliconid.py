"""
Pairs the FASTQ files and compresses it into a j3x file and j3x.stats file
j3x file output includes translocations and the mutations
j3x stats file is a summary of what is included in the j3x
"""

from manifestExtraction import grouper
from unpairedFastqProc import *
from ReadPairer import *
from ReadCompressor import *
from tqdm import tqdm
from multiprocessing import cpu_count
from statistics import median

from pprint import pprint

import os
import json
import time

class MainAmpliconID:
    def __init__(self):
        self.inDir = "data/1-raw/"
        self.outDir = "data/2-paired/"

        self.numThreads = cpu_count()
        self.chunksize = 250
        self.readPairer = ReadPairer()
        self.bytesPerRead = 350 # Estimated

    def test(self, inDir, outDir, configFile, referenceFile, filenameArray):
        self.inDir = inDir
        self.outDir = outDir
        self.readPairer = ReadPairer(configFile, referenceFile)
        
        self.processFiles(filenameArray)

    def run(self):
        print("MILo Amplicon Pairer")
        print("Chunksize (Process Pool): {0}".format(self.chunksize))
        print("Number of Threads: {0}".format(self.numThreads))
        print()
        
        with open('files.json') as file_list_file:    
            filenameArray = json.load(file_list_file)
            self.processFiles(filenameArray)

    def processFiles(self, filenameArray):
        for filenames in filenameArray:
            fq1, fq2, paired, skip = self.readFilenames(filenames)
            if skip:
                continue
            if (fq1 == '' or fq2 == '' or paired == ''):
                print('Please set the keys "fastq1", "fastq2" and "paired" in the files.json file')
                return
        
        for filenames in filenameArray:
            fq1, fq2, paired, skip = self.readFilenames(filenames)
            if skip:
                continue
            self.pairToJ3X(fq1, fq2, paired) 
        
    def readFilenames(self, filenames):
        fq1 = fq2 = paired = ''
        skip = False
        
        if ( 'fastq1' in filenames ):
            fq1 = filenames['fastq1']
        if ( 'fastq2' in filenames ):
            fq2 = filenames['fastq2']
        if ( 'paired' in filenames ):
            paired = filenames['paired']
        if ( 'skip' in filenames ):
            skip = filenames['skip']
            
        return fq1, fq2, paired, skip
        
    def pairToJ3X(self, fq1, fq2, paired):
        # Stats from read pairer
        totalReads = 0
        failedToPairReads = 0
        failedToMatchReads = 0
        failedToPairMatchReads = 0
        
        readCompressor = ReadCompressor(self.readPairer.getReferenceCount())
        if not os.path.exists(self.outDir):
            os.makedirs(self.outDir)
        with open(self.inDir + fq1) as fq1File, open(self.inDir + fq2) as fq2File:
            filesize1 = os.path.getsize(self.inDir + fq1)
            filesize2 = os.path.getsize(self.inDir + fq1)
            filesize = (filesize1 + filesize2) / 2
            estimatedReads = int(filesize / self.bytesPerRead)

            outfile = self.outDir + paired
            with open(outfile, "w+", newline = "") as outFile:
                print("{0} Crunching {1} and {2}".format(time.strftime('%X %d %b %Y'), fq1, fq2))
                start = time.time()
                
                # Creates iterators which deliver the 4 lines of each FASTQ read as a zip (ID, Sequence, Blank, Quality)
                fq1Iter, fq2Iter = grouper(fq1File, 4), grouper(fq2File, 4)
                processManager = ProcessPoolExecutor(self.numThreads)
                with tqdm(total=estimatedReads) as pbar:
                    result = processManager.map(self.readPairer.alignAndMerge, fq1Iter, fq2Iter, chunksize = self.chunksize)
                    for i, data in tqdm(enumerate(result)):
                        failedToPair = data.failedToPair
                        matchType = data.matchType
                        totalReads += 1
                        failedToPairReads += failedToPair
                        failedToMatchReads += 1 if matchType == 'nah' else 0
                        if failedToPair == 1 and matchType == 'nah':
                            failedToPairMatchReads += 1
                            
                        readCompressor.putPairedRead(data)
                        # outFile.write(data)
                        # outFile.write("\n\n")
                        pbar.update()
                pbar.close()
            
                readCompressor.prepareForCompression()
                print("\nCompressing\n")
                
                j3xSeqs = readCompressor.getDataList()
                
                # Print all the discarded stuff into another file
                discardedList = readCompressor.getDiscardedList()
                discardOutDir = self.outDir + 'discarded/'
                discardFile = discardOutDir + paired
                if not os.path.exists(discardOutDir):
                    os.makedirs(discardOutDir)
                with open(discardFile, "w+", newline = "") as discardedOutFile:
                    for discardSeq in discardedList:
                        self.writeSeqToFile(discardedOutFile, discardSeq)
                    discardedOutFile.close()
                
                numMergeAttempts, mergedCount, mergedUnsureCount, mergedD1Count, mergedD2Count, discardCountList, ampliconCountList, failedMergeAndDiscarded, templateCount, templateCountList = readCompressor.getStats()
                
                # Calculate the total number of discards, and the discard rates relative to each amplicon's read depth
                numDiscarded = sum(discardCountList)
                discardRates = [discarded / total if total != 0 else None for discarded, total in zip(discardCountList, ampliconCountList)]
                avgDiscardRate = median([rate for rate in discardRates if rate != None]) * 100

                # Format and write to j3x
                numUsableReads = 0
                for seq in j3xSeqs:
                    numReads = seq[1][0]
                    numUsableReads += numReads
                    self.writeSeqToFile(outFile, seq)

                outFile.close()
                
                # Calculate j3x statistics
                numUsableReads = int(numUsableReads)
                numSeqs = len(j3xSeqs)
                numOriginal = numDiscarded + numUsableReads
                prcntCompression = 100 - self.perc(numSeqs, numOriginal)
                prcntUsable = self.perc(numUsableReads, numOriginal)
                prcntMerged = self.perc(mergedCount , numOriginal)
                prcntMergedUnsure = self.perc(mergedUnsureCount , numOriginal)
                prcntDiscarded = self.perc(numDiscarded , numOriginal)
                prcntFailedPair = self.perc(failedToPairReads, numOriginal)
                prcntFailedMatch = self.perc(failedToMatchReads, numOriginal)
                prcntFailedPairMatch = self.perc(failedToPairMatchReads, numOriginal)
                
                with open(outfile + ".stats", "w+", newline = "") as statsFile:
                    statsFile.write("Overall File Stats\n")
                    pwrite(statsFile, "Specimen:        , {0}\t".format(paired))
                    pwrite(statsFile, "TimeTaken / Time:, {0}s\t, {1}\t".format(self.niceRound(time.time() - start), time.strftime('%X %d %b %Y')))
                    pwrite(statsFile, "Original Reads:  , {0}\t".format(numOriginal))
                    
                    pwrite(statsFile, "Pair pass:       , {0}\t, {1}%\t".format(numOriginal - failedToPairReads, 100 - prcntFailedPair))
                    pwrite(statsFile, "Pair fail:       , {0}\t, {1}%\t".format(failedToPairReads, prcntFailedPair))
                    pwrite(statsFile, "Match pass:      , {0}\t, {1}%\t".format(numOriginal - failedToMatchReads, 100 - prcntFailedMatch))
                    pwrite(statsFile, "Match fail:      , {0}\t, {1}%\t".format(failedToMatchReads, prcntFailedMatch))
                    pwrite(statsFile, "Pair/Match pass: , {0}\t, {1}%\t".format(numOriginal - failedToPairMatchReads, 100 - prcntFailedPairMatch))
                    pwrite(statsFile, "Pair/Match fail: , {0}\t, {1}%\t".format(failedToPairMatchReads, prcntFailedPairMatch))
                    
                    pwrite(statsFile, "Compressed Seq:  , {0}\t".format(numSeqs))
                    pwrite(statsFile, "Compression:     , {0}%\t".format(prcntCompression))
                    pwrite(statsFile, "Templates:       , {0}\t".format(templateCount))
                    pwrite(statsFile, "Usable Data:     , {0}\t, {1}%\t".format(numUsableReads, prcntUsable))
                    pwrite(statsFile, "Merged Data:     , {0}\t, {1}%\t, Data that was similar to templates".format(mergedCount, prcntMerged))
                    pwrite(statsFile, "Merged >1 Data:  , {0}\t, {1}%\t, Merges with more than one possible template candidate".format(mergedUnsureCount, prcntMergedUnsure))
                    pwrite(statsFile, "Merges D1:       , {0}\t, {1}%\t, Merges that had distance 1".format(mergedD1Count, self.perc(mergedD1Count, mergedCount)))
                    pwrite(statsFile, "Merges D2:       , {0}\t, {1}%\t, Merges that had distance 2".format(mergedD2Count, self.perc(mergedD2Count, mergedCount)))
                    pwrite(statsFile, "Discarded Data:  , {0}\t, {1}%\t".format(numDiscarded, prcntDiscarded))
                    pwrite(statsFile, "!Merge & Discard , {0}\t, , Discarded due to failed merge".format(failedMergeAndDiscarded))
                    pwrite(statsFile, "Avg Discard Rate:, {0}%\t, , Average discard rate per amplicon".format(self.niceRound(avgDiscardRate)))

                    statsFile.write("\nReference Amplicon Stats\n")
                    statsFile.write("{0}, {1}, {2}, {3}, {4}\n".format('Ref ampID', 'Amplicon Count', 'Template Count', 'Discard Count', 'Discard Rate'))
                    for i in range(len(templateCountList)):
                        discardPercent = '-'
                        if discardRates[i] != None:
                            discardPercent = self.niceRound(discardRates[i])
                        statsFile.write("{0}, {1}, {2}, {3}, {4}\n".format(i, ampliconCountList[i], templateCountList[i], discardCountList[i], discardPercent))
                    
                    # statsFile.close()
                    
    def writeSeqToFile(self, oFile, seq):
        sequence = seq[0]
        numReads = seq[1][0]
        numReadsMerged = seq[1][1]
        infoLine = seq[1][2]
        quality = seq[1][3]
        oFile.write("{0}, R:{1}, M:{2}".format(infoLine, int(numReads), int(numReadsMerged)))
        oFile.write('\n')
        oFile.write(sequence)
        oFile.write('\n')
        oFile.write(quality)
        oFile.write("\n\n")
        
    def niceRound(self, num):
        return int(num * 100)/100
    def perc(self, numerator, denominator):
        if denominator == 0: return 0
        return int(numerator * 1000 / denominator) / 10

if __name__ ==  "__main__":
    main = MainAmpliconID()
    main.run()