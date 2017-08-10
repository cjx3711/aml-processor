import sys
from j4xUtils import *

if ( len(sys.argv) < 3 ):
    print ("To use this script, please provide 2 arguments. Base sequence and mutate hash")
    print ('E.g.')
    print ('> python sequenceMutator.py "AAATTTCCCGGG" "S:0:A-T I:9:-TT D:11:G-"')
    print ('> TAATTTCCCTTGG')

else:

    base = sys.argv[1]
    mutation = sys.argv[2]
    mutated = mutateSequenceByHash(base, mutation)

    print(mutated)