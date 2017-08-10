import glob
import os
from j4xUtils import *
from pprint import pprint

inDir = "data/3-mutations/"
outDir = "data/4-mutationstats"

filenaneEnd = "_MUTATIONS.j4x"
mutationHumanDictionary = {}

def run():
    minOccurences = 2
    for filepath in glob.glob(os.path.join(inDir, '*.j4x')):
        # Read all files in input directory
        filenameParts = filepath.split("/")
        filenameOnly = filenameParts[len(filenameParts)-1]
        if (filenameOnly[-len(filenaneEnd):] == filenaneEnd):
            personName = filenameOnly[:-len(filenaneEnd)]
            proessSingleHuman(filepath, personName)
    
    mutationTupleList = list(mutationHumanDictionary.items())
    filteredTupleList = [x for x in mutationTupleList if x[1]['fileOccurrences'] >= minOccurences]
    
    filteredTupleList.sort(key=lambda tup: -tup[1]['totalOccurrences'])
    pprint(filteredTupleList)
    
        

def proessSingleHuman(filepath, personName):
    with open(filepath) as mutationFile:
        for line in mutationFile:
            processSingleLine(line[:-1], personName) # Remove the \n

def processSingleLine(line, personName):
    parts = line.split(", ")
    occurences = int(parts[0])
    mutationHash = parts[1]
    
    if ( mutationHash not in mutationHumanDictionary ):
        mutationHumanDictionary[mutationHash] = {
            'totalOccurrences' : 0,
            'fileOccurrences': 0,
            'occurrences': [],
            'humans': []
        }
    
    mutationHumanDictionary[mutationHash]['totalOccurrences'] += occurences
    mutationHumanDictionary[mutationHash]['fileOccurrences'] += 1
    mutationHumanDictionary[mutationHash]['occurrences'].append(occurences)
    mutationHumanDictionary[mutationHash]['humans'].append(personName)
    
run()