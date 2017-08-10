from manifestExtraction import grouper
from unpairedFastqProc import *
from ReadPairer import *
from tqdm import tqdm
from multiprocessing import cpu_count

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
    with open('files.json') as file_list_file:    
        filenameArray = json.load(file_list_file)
    
        for filenames in filenameArray:
            fq1, fq2, paired = readFilenames(filenames)
            
            if (fq1 == '' or fq2 == '' or paired == ''):
                print('Please set the keys "fastq1", "fastq2" and "paired" in the config.json file')
                return
        
        for filenames in filenameArray:
            fq1, fq2, paired = readFilenames(filenames)
    
            # Don't understand the __name__ thing, but it's required according to SO
            if __name__ ==  "__main__": pairToJ3X(fq1, fq2, paired, inDir, outDir) 
    
    
def readFilenames(filenames):
    fq1 = fq2 = paired = ''
    
    if ( 'fastq1' in filenames ):
        fq1 = filenames['fastq1']
    if ( 'fastq2' in filenames ):
        fq2 = filenames['fastq2']
    if ( 'paired' in filenames ):
        paired = filenames['paired']
    return fq1, fq2, paired
    
def pairToJ3X(fq1, fq2, paired, inDir, outDir):
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
                    outFile.write(data)
                    outFile.write("\n\n")
                    pbar.update()
                outFile.close()
            pbar.close()

            print("{0} Dumped {1}".format(time.strftime('%X %d %b %Y'), paired))
            print("Took {0}s".format(time.time() - start))

run()