from difflib import SequenceMatcher

def mutationID(base, compare):
    seqs = SequenceMatcher(None, base, compare, autojunk = False)
    matches = seqs.get_matching_blocks();
    start = [0, 0]
    mutations = []
    
    for i in range(0, len(matches)):
        
        # Get the bits from the start to the current point
        posA = matches[i].a
        posB = matches[i].b
        deltaA = posA - start[0]
        deltaB = posB - start[1]
        
        if ( deltaA == 0 and deltaB == 0 ):
            # No changes
            pass
        elif ( deltaA == deltaB ):
            # Substitution
            # print('Substitution at {0} from {1} to {2}'.format(posA-deltaA, base[posA-deltaA:posA], compare[posB-deltaB:posB]))
            mutations.append({ 'type': 'S', 'pos': posA-deltaA, 'from': base[posA-deltaA:posA], 'to': compare[posB-deltaB:posB]})
        elif ( deltaA > deltaB ):
            if ( deltaB == 0 ):
                # Simple deletion
                # print('Deletion at {0} of {1}'.format(posA-deltaA, base[posA-deltaA:posA]))
                mutations.append({ 'type': 'D', 'pos': posA-deltaA, 'from': base[posA-deltaA:posA], 'to': ''})
                
            else:
                # Deletion and insertion
                # print('DelAdd at {0} del {1} add {2}'.format(posA-deltaA, base[posA-deltaA:posA], compare[posB-deltaB:posB]))
                mutations.append({ 'type': 'S', 'pos': posA-deltaA, 'from': base[posA-deltaA:posA], 'to': compare[posB-deltaB:posB]})
                
        elif ( deltaA < deltaB ):
            # Addition
            if ( deltaA == 0 ):
                # Simple Addition
                # print('Addition at {0} of {1}'.format(posA-deltaA, compare[posB-deltaB:posB])) 
                mutations.append({ 'type': 'I', 'pos': posA-deltaA, 'from': '', 'to': compare[posB-deltaB:posB]})
                
            else:
                # Insertion and Deletion
                # print('AddDel at {0} add {1} del {2}'.format(posA-deltaA, compare[posB-deltaB:posB], base[posA-deltaA:posA]))
                mutations.append({ 'type': 'S', 'pos': posA-deltaA, 'from': base[posA-deltaA:posA], 'to': compare[posB-deltaB:posB]})
                
                
        
        start[0] = posA + matches[i].size
        start[1] = posB + matches[i].size
    # print(mutations)
    # print(mutationArrayToHash(mutations))
    
    return mutations
    
    # print(matches[i])
    # print(deltaA)
    # print(deltaB)
def mutationArrayToHash(mutations):
    mutationStrings = []
    for m in mutations:
        mutationStrings.append(mutationToString(m))
    return ' '.join(mutationStrings)

def mutationIDAsHash(base, compare):
    return mutationArrayToHash(mutationID(base, compare))
    
def mutationToString(mutation):
    return mutation['type'] + ':' + str(mutation['pos']) + ':' + str(mutation['from']) + '-' + str(mutation['to'])
    
def hashToMutationArray(mutationHash):
    mutations = mutationHash.split(' ')
    mutationArray = []
    for mutation in mutations:
        if (mutation.startswith(('S', 'I', 'D'))):
            mutationArray.append(hashToMutation(mutation))
    
    return mutationArray
    
def hashToMutation(mutationHash):
    hashParts = mutationHash.split(':')
    mutationType = hashParts[0]
    mutationPosition = int(hashParts[1])
    changeParts = hashParts[2].split('-')
    mutationFrom = changeParts[0]
    mutationTo = changeParts[1]
    return {
        'type': mutationType,
        'pos': mutationPosition,
        'from': mutationFrom,
        'to': mutationTo
    }

def convertHashPositionsToCoordinates(mutationHash, startCoordinate):
    positionArray = extractPositionsFromHash(mutationHash)
    coordinateArray = [ str(x + startCoordinate) for x in positionArray ]
    return ' '.join(coordinateArray)
    
def extractPositionsFromHash(mutationHash):
    mutationArray = hashToMutationArray(mutationHash)
    positionArray = []
    for mutation in mutationArray:
        positionArray.append(mutation['pos'])
    return positionArray
    
def mutateSequenceByHash(base, mutationHash):
    mutationArray = hashToMutationArray(mutationHash)
    return mutateSequence(base, mutationArray)
    
def mutateSequence(base, mutationArray):
    cursor = 0
    mutated = ''
    for mutation in mutationArray:
        mutated += base[cursor:mutation['pos']]
        cursor = mutation['pos']
        if ( mutation['type'] == 'I' ):
            mutated += mutation['to']
        elif ( mutation['type'] == 'D' ):
            actual = base[cursor:cursor+len(mutation['from'])]
            if ( actual != mutation['from'] ):
                return 'Invalid mutation'
            cursor += len(mutation['from'])
            
        elif ( mutation['type'] == 'S' ):
            actual = base[cursor:cursor+len(mutation['from'])]
            if ( actual != mutation['from'] ):
                return 'Invalid mutation'
            cursor += len(mutation['from'])
            mutated += mutation['to']
            
    mutated += base[cursor:]
    return mutated