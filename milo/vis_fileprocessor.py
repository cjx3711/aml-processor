"""
Modelled after main_ampliconid.py
Pairs the FASTQ files and runs calculations on the reads for visualisations

"""

from manifestExtraction import grouper
from unpairedFastqProc import *
from ReadPairAndID import *
from ReadCompressor import *
from ProbReadCompressor import *
from tqdm import tqdm
from multiprocessing import cpu_count
from genomicsUtils import reverseComplement
from generalUtils import *
from MergeVisualisationMetrics import *
from pprint import pprint

import os
import json
import time
import csv
from collections import defaultdict
from statistics import median


class VisFileProcessor:
    def __init__(self):
        self.numThreads = cpu_count()
        self.chunksize = 250
        self.bytesPerRead = 350  # Estimated
        self.configFile = 'config.json'
        self.referenceFile = 'references/Manifest.csv'
        with open(self.referenceFile) as refFile:
            self.directionList = ["+"] + [line[2] for line in csv.reader(refFile, delimiter=',')][1:]
        self.conflictingDirTrans = defaultdict(int)

    def test(self, inDir, outDir, configFile, referenceFile, filenameArray):
        self.configFile = configFile
        self.inDir = inDir
        self.outDir = outDir
        self.readPairer = ReadPairAndID(configFile, referenceFile)

        for i, filenames in enumerate(filenameArray):
            print("\n=================================================")
            print("Processing {0}. {1}/{2}".format(filenames.fastq1, i, len(filenameArray)))
            self.pairToJ3X(filenames.fastq1, filenames.fastq2, filenames.paired, filenames.pairedStats)

    def run(self):
        self.inDir = "data/1-raw/"
        self.outDir = "data/2-paired/"

        self.readPairer = ReadPairAndID()

        print("MILo Amplicon Pairer")
        print("Chunksize (Process Pool): {0}".format(self.chunksize))
        print("Number of Threads: {0}".format(self.numThreads))
        print()

        filenameArray = getFileList('files.json')
        for i, filenames in enumerate(filenameArray):
            print("\n=================================================")
            print("Processing {0}. {1}/{2}".format(filenames.fastq1, i, len(filenameArray)))
            self.pairToJ3X(filenames.fastq1, filenames.fastq2, filenames.paired, filenames.pairedStats)

    def pairToJ3X(self, fq1, fq2, paired, pairedStats):
        # Stats from read pairer
        totalReads = 0

        mergeVisualisation = MergeVisualisationMetrics()
        if not os.path.exists(self.outDir):
            os.makedirs(self.outDir)
        with open(self.inDir + fq1) as fq1File, open(self.inDir + fq2) as fq2File:
            filesize1 = os.path.getsize(self.inDir + fq1)
            filesize2 = os.path.getsize(self.inDir + fq1)
            filesize = (filesize1 + filesize2) / 2
            estimatedReads = int(filesize / self.bytesPerRead)

            outfile = self.outDir + paired
            with open(outfile, "w+", newline="") as outFile:
                print("{0} Crunching {1} and {2}".format(time.strftime('%X %d %b %Y'), fq1, fq2))
                start = time.time()

                # Creates iterators which deliver the 4 lines of each FASTQ read as a zip (ID, Sequence, Blank, Quality)
                fq1Iter, fq2Iter = grouper(fq1File, 4), grouper(fq2File, 4)
                processManager = ProcessPoolExecutor(self.numThreads)
                with tqdm(total=estimatedReads) as pbar:
                    # Align unpaired reads, merge them, and identify their amplicon
                    result = processManager.map(self.readPairer.alignAndMerge, fq1Iter, fq2Iter,
                                                chunksize=self.chunksize)
                    for i, pairedRead in tqdm(enumerate(result)):
                        totalReads += 1

                        mergeVisualisation.putPairedRead(pairedRead)
                        pbar.update()
                pbar.close()

                mergeVisualisation.closeFile()
                print("Done")


if __name__ == "__main__":
    main = VisFileProcessor()
    main.run()