from manifestExtraction import grouper
from unpairedFastqProc import *
from ReadPairer import *
from ReadCompressor import *
from tqdm import tqdm
from multiprocessing import cpu_count

from pprint import pprint

import os
import json
import time


inDir = "data/1-raw/"
outDir = "data/2-paired/"

numThreads = cpu_count()
chunksize = 250
readPairer = ReadPairer(probabilistic = False)
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
            # pool = Pool(numThreads)
            processManager = ProcessPoolExecutor(numThreads)
            with tqdm(total=estimatedReads) as pbar:
                result = processManager.map(readPairer.alignAndMerge, fq1Iter, fq2Iter, chunksize = chunksize)
                for i, data in tqdm(enumerate(result)):
                    readCompressor.putPairedRead(data)
                    # outFile.write(data)
                    # outFile.write("\n\n")
                    pbar.update()
            pbar.close()
            
            
            
            sortedCompressedList, totalOnes, matchedOnes = readCompressor.getDataList()
            discardedOnes = totalOnes - matchedOnes
            totalMatched = 0
            for read in sortedCompressedList:
                sequence = read[0]
                count = read[1][0]
                matchCount = read[1][1]
                iddata = read[1][2]
                quality = read[1][3]
                totalMatched += count
                outFile.write("{0}, C:{1}, M:{2}".format(iddata, count, matchCount))
                outFile.write('\n')
                outFile.write(sequence)
                outFile.write('\n')
                outFile.write(quality)
                outFile.write("\n\n")
            outFile.close()
            
            lines = len(sortedCompressedList)
            total = discardedOnes + totalMatched
            percCompression = 100 - int(lines / total * 1000)/10
            percUsable = int(totalMatched / total * 1000)/10
            percMatched = int(matchedOnes / total * 1000)/10
            percDiscarded = int(discardedOnes / total * 1000)/10
            
            with open(outfile + ".stats", "w+", newline = "") as statsFile:
                pwrite(statsFile, "{0} Dumped {1}".format(time.strftime('%X %d %b %Y'), paired))
                pwrite(statsFile, "Compressed {0} reads into {1} lines. ({2}%% compression)".format(total, lines, percCompression))
                pwrite(statsFile, "Usable data on {0} of {1}. ({2}%)".format(totalMatched, total, percUsable))
                pwrite(statsFile, "Close match on {0} of {1}. ({2}%)".format(matchedOnes, total, percMatched))
                pwrite(statsFile, "Discarded {0} of {1}. ({2}%)".format(discardedOnes, total, percDiscarded))
                pwrite(statsFile, "Took {0}s\n\n".format(time.time() - start))
                
def pwrite(file, message):
    print(message)
    file.write(message + "\n")
if __name__ ==  "__main__": run()