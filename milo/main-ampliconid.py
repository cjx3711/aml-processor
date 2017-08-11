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
            
            sortedCompressedList = readCompressor.getDataList()
            for read in sortedCompressedList:
                sequence = read[0]
                count = read[1][0]
                iddata = read[1][1]
                quality = read[1][2]
                outFile.write("{0}, {1}".format(iddata, count))
                outFile.write('\n')
                outFile.write(sequence)
                outFile.write('\n')
                outFile.write(quality)
                outFile.write("\n\n")
            outFile.close()
                
            print("{0} Dumped {1}".format(time.strftime('%X %d %b %Y'), paired))
            print("Took {0}s\n\n".format(time.time() - start))

if __name__ ==  "__main__": run()