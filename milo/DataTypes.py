from collections import namedtuple

# Used for file reading
FileTypes = namedtuple('FileTypes', ['fastq1', 'fastq2', 'paired', 'pairedStats', 'tiled', 'mutation'])


# Used for TranslocatedBlockMatcher
EmptyMatch = namedtuple('EmptyMatch', ['a', 'b', 'size'])

# Used for the return value of ReadPairer.alignAndMerge
# DEPRECATED. Replaced by PairedRead
AlignedAndMerged = namedtuple('AlignedAndMerged', ['failedToPair', 'matchType', 'sequenceInfo', 'sequenceData', 'qualityData'])

# Used for return value of ReadPairAndID.alignAndMerge
PairedRead = namedtuple('PairedRead', ['baseSeq', 'phredSeq', 'pairSuccess'])

# Used to hold the j3x reference matrix
J3xStats = namedtuple('J3xStats', ['ampID', 'totalReads', 'numTemplates', 'numDiscards'])

# Used by the GeneMap
GeneReference = namedtuple('GeneReference', ['geneID', 'count', 'ampIDs', 'name', 'sequence', 'firstCoord', 'lastCoord', 'startCoords', 'endCoords'])