import glob
import os
from j4xUtils import *
from pprint import pprint

inDir = os.path.join("data", "3-mutations")
outDir = os.path.join("data", "4-mutationstats")
significanceThreshold = 0.2

filenameEnd = "_MUTATIONS.j4x"
mutationHumanDictionary = {}

def run():
    minOccurences = 1
    for filepath in glob.glob(os.path.join(inDir, '*.j4x')):
        # Read all files in input directory
        filenameParts = filepath.split(os.sep)
        filenameOnly = filenameParts[len(filenameParts)-1]
        if (filenameOnly[-len(filenameEnd):] == filenameEnd):
            personName = filenameOnly[:-len(filenameEnd)]
            processSingleHuman(filepath, personName)
    
    mutationTupleList = list(mutationHumanDictionary.items())
    filteredTupleList = [x for x in mutationTupleList if x[1]['fileOccurrences'] >= minOccurences]
    
    # Stably sort the list by file occurrences,
    filteredTupleList.sort(key = lambda tup: tup[1]['fileOccurrences'], reverse = True)
    # Average VAFrequency
    filteredTupleList.sort(key = lambda tup: sum(tup[1]['VAFrequency']) / tup[1]['fileOccurrences'], reverse = True)
    # 
    filteredTupleList.sort(key = lambda tup: tup[0][:tup[0].find(" ")])
    significantTupleList = [("",)]
    for x in filteredTupleList:
        if ( sum(x[1]['VAFrequency']) / x[1]['fileOccurrences'] > significanceThreshold ):
            significantTupleList.append(x)
        elif ( x[0][:x[0].find(" ")] == significantTupleList[-1][0][:significantTupleList[-1][0].find(" ")]):
            significantTupleList.append(x)
    pprint(significantTupleList)
    
        

def processSingleHuman(filepath, personName):
    with open(filepath) as mutationFile:
        for line in mutationFile:
            if line[:-1] == "Reference Amplicon Stats":
                break
            processSingleLine(line[:-1], personName) # Remove the \n

def processSingleLine(line, personName):
    parts = line.split(", ")
    readSplitPt = parts[0].find("/")
    occurrences = int(parts[0][:readSplitPt])
    vaFrequency = round(occurrences / int(parts[0][readSplitPt + 1:]), 2)
    mutationHash = parts[1]
    
    if ( mutationHash not in mutationHumanDictionary ):
        mutationHumanDictionary[mutationHash] = {
            'totalOccurrences' : 0,
            'fileOccurrences': 0,
            'occurrences': [],
            'VAFrequency': [],
            'humans': []
        }
    
    mutationHumanDictionary[mutationHash]['totalOccurrences'] += occurrences
    mutationHumanDictionary[mutationHash]['fileOccurrences'] += 1
    mutationHumanDictionary[mutationHash]['occurrences'].append(occurrences)
    mutationHumanDictionary[mutationHash]['VAFrequency'].append(vaFrequency)
    mutationHumanDictionary[mutationHash]['humans'].append(personName)
    
run()