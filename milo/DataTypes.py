from collections import namedtuple

# Used for TranslocatedBlockMatcher
EmptyMatch = namedtuple('EmptyMatch',
                        ['a', 'b', 'size']
                       )

# Used for the return value of ReadPairer.alignAndMerge
AlignedAndMerged = namedtuple('AlignedAndMerged',
                             ['failedToPair', 'matchType', 'sequenceInfo', 'sequenceData', 'qualityData' ]
                             )