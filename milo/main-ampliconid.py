from manifestExtraction import *
from unpairedFastqProc import *
from ReadPairer import *
import json

inDir = "data/Raw/"
outDir = "data/Processed/"

numThreads = 12
readPairer = ReadPairer(probabilistic = False)

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
        outfile = outDir + paired
        print("Crunching {0} and {1}".format(fq1, fq2))
        
        with open(outfile, "w+", newline = "") as outFile:
            # Creates iterators which deliver the 4 lines of each FASTQ read as a zip (ID, Sequence, Blank, Quality)
            fq1Iter, fq2Iter = grouper(fq1File, 4), grouper(fq2File, 4)
            with ProcessPoolExecutor(numThreads) as processManager:
                # Calls alignAndMerge(FASTQ1's (ID, Sequence, Blank, Quality), FASTQ2's (ID, Sequence, Blank, Quality))
                for x in processManager.map(readPairer.alignAndMerge, fq1Iter, fq2Iter, chunksize = 250):
                    outFile.write(x)
                    outFile.write("\n\n")
                outFile.close()
                print("Dumped {0}".format(paired))

run()