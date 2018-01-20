from difflib import SequenceMatcher
import re

sequenceMatcher = SequenceMatcher(None, '', '', autojunk = False)
ref = ''

def mutationID(refInput, read):
    global ref
    if ( ref != refInput ):
        ref = refInput;
        sequenceMatcher.set_seq2(ref)
    sequenceMatcher.set_seq1(read)
    matches = sequenceMatcher.get_matching_blocks()
    start = [0, 0]
    mutations = []
    
    for i in range(0, len(matches)):
        # Get the bits from the start to the current point
        posRef = matches[i].b
        posRead = matches[i].a
        deltaRef = posRef - start[0]
        deltaRead = posRead - start[1]
        
        if ( deltaRef == 0 and deltaRead == 0 ):
            # No changes
            pass
        elif ( deltaRef == deltaRead ):
            # Substitution
            # print('Substitution at {0} from {1} to {2}'.format(posRef-deltaRef, refInput[posRef-deltaRef:posRef], read[posRead-deltaRead:posRead]))
            mutations.append({ 'type': 'S', 'pos': posRef-deltaRef, 'from': refInput[posRef-deltaRef:posRef], 'to': read[posRead-deltaRead:posRead]})
        elif ( deltaRef > deltaRead ):
            if ( deltaRead == 0 ):
                # Simple deletion
                # print('Deletion at {0} of {1}'.format(posRef-deltaRef, refInput[posRef-deltaRef:posRef]))
                mutations.append({ 'type': 'D', 'pos': posRef-deltaRef, 'from': refInput[posRef-deltaRef:posRef], 'to': ''})
                
            else:
                # Deletion and insertion
                # print('DelAdd at {0} del {1} add {2}'.format(posRef-deltaRef, refInput[posRef-deltaRef:posRef], read[posRead-deltaRead:posRead]))
                mutations.append({ 'type': 'S', 'pos': posRef-deltaRef, 'from': refInput[posRef-deltaRef:posRef], 'to': read[posRead-deltaRead:posRead]})
                
        elif ( deltaRef < deltaRead ):
            # Addition
            if ( deltaRef == 0 ):
                # Simple Addition
                # print('Addition at {0} of {1}'.format(posRef-deltaRef, read[posRead-deltaRead:posRead])) 
                mutations.append({ 'type': 'I', 'pos': posRef-deltaRef, 'from': '', 'to': read[posRead-deltaRead:posRead]})
                
            else:
                # Insertion and Deletion
                # print('AddDel at {0} add {1} del {2}'.format(posRef-deltaRef, read[posRead-deltaRead:posRead], refInput[posRef-deltaRef:posRef]))
                mutations.append({ 'type': 'S', 'pos': posRef-deltaRef, 'from': refInput[posRef-deltaRef:posRef], 'to': read[posRead-deltaRead:posRead]})
                
                
        
        start[0] = posRef + matches[i].size
        start[1] = posRead + matches[i].size
    # print(mutations)
    # print(mutationArrayToHash(mutations))
    
    return mutations
    
    # print(matches[i])
    # print(deltaRef)
    # print(deltaRead)
def mutationArrayToHash(mutations):
    mutationStrings = []
    for m in mutations:
        mutationStrings.append(mutationToString(m))
    return ' '.join(mutationStrings)

def mutationIDAsHash(refInput, read):
    return mutationArrayToHash(mutationID(refInput, read))
    
def mutationToString(mutation):
    return mutation['type'] + ':' + str(mutation['pos']) + ':' + str(mutation['from']) + '-' + str(mutation['to'])
    
def hashToMutationArray(mutationHash):
    mutationList = re.split("([SID])", mutationHash)
    mutationList = [mutantType + _ for mutantType, _ in zip(mutationList[1::2], mutationList[::2][1:])]
    return [hashToMutation(x) for x in ([mutationList[i][:-1] for i in range(len(mutationList) - 1)] + [mutationList[-1]])]
    
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
    
def mutateSequenceByHash(refInput, mutationHash):
    mutationArray = hashToMutationArray(mutationHash)
    return mutateSequence(refInput, mutationArray)
    
def mutateSequence(refInput, mutationArray):
    cursor = 0
    mutated = ''
    for mutation in mutationArray:
        mutated += refInput[cursor:mutation['pos']]
        cursor = mutation['pos']
        if ( mutation['type'] == 'I' ):
            mutated += mutation['to']
        elif ( mutation['type'] == 'D' ):
            actual = refInput[cursor:cursor+len(mutation['from'])]
            if ( actual != mutation['from'] ):
                return 'Invalid mutation'
            cursor += len(mutation['from'])
            
        elif ( mutation['type'] == 'S' ):
            actual = refInput[cursor:cursor+len(mutation['from'])]
            if ( actual != mutation['from'] ):
                return 'Invalid mutation'
            cursor += len(mutation['from'])
            mutated += mutation['to']
            
    mutated += refInput[cursor:]
    return mutated