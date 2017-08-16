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


inDir = "data/1-raw/"
outDir = "data/2-paired/"

numThreads = cpu_count()
chunksize = 250
readPairer = ReadPairer()
bytesPerRead = 350 # Estimated

def run():
    print("MILo Amplicon Pairer")
    print("Chunksize (Process Pool): {0}".format(chunksize))
    print()
    
    with open('files.json') as file_list_file:    
        filenameArray = json.load(file_list_file)
    
        for filenames in filenameArray:
            fq1, fq2, paired, skip = readFilenames(filenames)
            if skip:
                continue
            if (fq1 == '' or fq2 == '' or paired == ''):
                print('Please set the keys "fastq1", "fastq2" and "paired" in the config.json file')
                return
        
        for filenames in filenameArray:
            fq1, fq2, paired, skip = readFilenames(filenames)
            if skip:
                continue
            # Don't understand the __name__ thing, but it's required according to SO
            pairToJ3X(fq1, fq2, paired, inDir, outDir) 

def readFilenames(filenames):
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
    
def pairToJ3X(fq1, fq2, paired, inDir, outDir):
    
    readCompressor = ReadCompressor()
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    with open(inDir + fq1) as fq1File, open(inDir + fq2) as fq2File:
        filesize1 = os.path.getsize(inDir + fq1)
        filesize2 = os.path.getsize(inDir + fq1)
        filesize = (filesize1 + filesize2) / 2
        estimatedReads = int(filesize / bytesPerRead)

        outfile = outDir + paired
        with open(outfile, "w+", newline = "") as outFile:
            print("{0} Crunching {1} and {2}".format(time.strftime('%X %d %b %Y'), fq1, fq2))
            start = time.time()
            
            # Creates iterators which deliver the 4 lines of each FASTQ read as a zip (ID, Sequence, Blank, Quality)
            fq1Iter, fq2Iter = grouper(fq1File, 4), grouper(fq2File, 4)
            processManager = ProcessPoolExecutor(numThreads)
            with tqdm(total=estimatedReads) as pbar:
                result = processManager.map(readPairer.alignAndMerge, fq1Iter, fq2Iter, chunksize = chunksize)
                for i, data in tqdm(enumerate(result)):
                    readCompressor.putPairedRead(data)
                    # outFile.write(data)
                    # outFile.write("\n\n")
                    pbar.update()
            pbar.close()
            
            # Counts the read depth of each amplicon
            ampliconCounts = [0] * readPairer.getReferenceCount()
            rawDataList = readCompressor.getRawDataList()
            for rawTuple in rawDataList:
                ampID = int(rawTuple[1][2].split(',')[0].strip()[3:]) # Extracts the ampID from the info line
                count = rawTuple[1][0]
                ampliconCounts[ampID] += count
            
            j3xSeqs, numMergeAttempts, mergedCount, mergedUnsureCount, mergedD1, mergedD2, discardCountList = readCompressor.getDataList(ampliconCounts)
            # Calculate the total number of discards, and the rates relative to each amplicon's read depth
            numDiscarded = sum(discardCountList)
            discardRates = [discarded / total if total != 0 else None for discarded, total in zip(discardCountList, ampliconCounts)]
            avgDiscardRate = median([x for x in discardRates if x != None])
            if any([abs(avgDiscardRate - x) > 0.1 for x in discardRates if x != None]): # If the discards are concentrated in one amplicon
                print("ALERT: The following amplicons have high discard rates:")
                print([(ampID, round(rate, 2)) for ampID, rate in enumerate(discardRates) if rate != None and abs(avgDiscardRate - rate) > 0.1])
                print("The average discard rate is {0}%".format(round(avgDiscardRate, 3) * 100))

            # Format and write to j3x
            totalAcrossAmplicons = 0
            for seq in j3xSeqs:
                sequence = seq[0]
                numReads = seq[1][0]
                numReadsMerged = seq[1][1]
                infoLine = seq[1][2]
                quality = seq[1][3]
                totalAcrossAmplicons += numReads
                outFile.write("{0}, R:{1}, M:{2}".format(infoLine, numReads, numReadsMerged))
                outFile.write('\n')
                outFile.write(sequence)
                outFile.write('\n')
                outFile.write(quality)
                outFile.write("\n\n")
            outFile.close()
            
            # Calculate j3x statistics
            numSeqs = len(j3xSeqs)
            numOriginal = numDiscarded + totalAcrossAmplicons
            prcntCompression = round(1 - numSeqs / numOriginal, 3) * 100
            prcntUsable = round(totalAcrossAmplicons / numOriginal, 3) * 100
            prcntMerged = round(mergedCount / numOriginal, 3) * 100
            prcntMergedUnsure = round(mergedUnsureCount / numOriginal, 3) * 100
            prcntDiscarded = round(numDiscarded / numOriginal, 3) * 100
            
            with open(outfile + ".stats", "w+", newline = "") as statsFile:
                pwrite(statsFile, "{0} Finished writing to j3x {1}".format(time.strftime('%X %d %b %Y'), paired))
                pwrite(statsFile, "Compressed {0} reads into {1} sequences. ({2}%% compression)".format(numOriginal, numSeqs, prcntCompression))
                pwrite(statsFile, "Usable data on {0} of {1}. ({2}%)".format(totalAcrossAmplicons, numOriginal, prcntUsable))
                pwrite(statsFile, "Merged {0} of {1}. ({2}%)".format(mergedCount, numOriginal, prcntMerged))
                pwrite(statsFile, "Match more than one on {0} of {1}. ({2}%)".format(mergedUnsureCount, numOriginal, prcntMergedUnsure))
                pwrite(statsFile, "Discarded {0} of {1}. ({2}%)".format(numDiscarded, numOriginal, prcntDiscarded))
                pwrite(statsFile, "{0} reads were merged with distance of 1, while {1} were merged with distance of 2.".format(mergedD1, mergedD2))
                pwrite(statsFile, "Took {0}s\n\n".format(time.time() - start))
                
def pwrite(file, message):
    print(message)
    file.write(message + "\n")
if __name__ ==  "__main__": run()