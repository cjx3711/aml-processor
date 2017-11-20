import sys
from j4xUtils import *
from genomicsUtils import *

if ( len(sys.argv) < 3 ):
    print ("To use this script, please provide 2 - 3 arguments. Base sequence and mutate hash")
    print ('E.g.')
    print ('> python sequenceMutator.py "AAATTTCCCGGG" "S:0:A-T I:9:-TT D:11:G-"')
    print ('> TAATTTCCCTTGG')
    print ('Put true as the third argument if you need to reverse complement the input')

else:

    base = sys.argv[1]
    mutation = sys.argv[2]
    if (len(sys.argv) == 4 and sys.argv[3] == 'true'):
        base = reverseComplement(base)
    mutated = mutateSequenceByHash(base, mutation)

    print(mutated)