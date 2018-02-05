"""
This is a helper class to process the reads and output some metrics on each read.
These metrics will be used for visualisation on the clustering of the reads.
"""

import csv
from collections import Counter
from phredUtils import *

class MergeVisualisationMetrics:

    def __init__(self, outputFile = 'vis.csv'):
        # Create file writer
        self.csvfile = open(outputFile, 'w')
        self.fileOut = csv.writer(self.csvfile)
        self.fileOut.writerow(['A', 'T', 'C', 'G', 'Variance', 'Mean', '20% Variance', '20% Mean', 'Delta Variance', 'Delta Mean'])
        
    def putPairedRead(self, pairedRead):
        if not pairedRead.pairSuccess:
            return

        baseSeq = pairedRead.baseSeq
        phredSeq = pairedRead.phredSeq
        phrednp = phredSeqToNpArray(phredSeq)
        # Process metrics
        # Get variance
        variance = np.var(phrednp)
        # Get average
        mean = np.mean(phrednp)
        # Get A, T, C, G count
        counter = Counter(baseSeq)
        aCount = counter['A']
        tCount = counter['T']
        cCount = counter['C']
        gCount = counter['G']

        # Get the bottom 20% of values
        k = max(int(phrednp.shape[0] * 0.2), 20)
        indices = np.argpartition(phrednp, k)
        lowerarray = phrednp[indices[:k]]

        lowvariance = np.var(lowerarray)
        lowmean = np.mean(lowerarray)

        deltavar = variance - lowvariance
        deltamean = mean - lowmean

        self.fileOut.writerow([aCount, tCount, cCount, gCount, variance, mean, lowvariance, lowmean, deltavar, deltamean])
        
    def closeFile(self):
        self.csvfile.close()