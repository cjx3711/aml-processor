import glob
import os
from j4xUtils import *
from pprint import pprint

statsDir = os.path.join("data", "2-paired")
inDir = os.path.join("data", "3-mutations")
outDir = os.path.join("data", "4-mutationstats")
significanceThreshold = 0.2

filenameEnd = "_MUTATIONS.j4x"
mutationHumanDictionary = {}


# Stats from the stats file of the reference amplicon
# Format:
# [
# (ampID, Amplicon Count, Template Count, Discard Count, Discard Rate),
# ... x 571
# ]
referenceAmpliconStats = [] 

def run():
    minOccurences = 1
    for filepath in glob.glob(os.path.join(inDir, '*.j4x')):
        # Read all files in input directory
        filenameParts = filepath.split(os.sep)
        filenameOnly = filenameParts[len(filenameParts)-1]
        if (filenameOnly[-len(filenameEnd):] == filenameEnd):
            personName = filenameOnly[:-len(filenameEnd)]
            statsFile = statsDir + os.sep + personName + '_PAIRED.j3x.stats' 
            processSingleHuman(filepath, statsFile, personName)
    
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
    
        

def processSingleHuman(filepath, statsfilepath, personName):
    print(statsfilepath)
    with open(statsfilepath) as statsFile:
        reading = False
        for line in statsFile:
            if reading == False:
                if line.startswith("Reference Amplicon Stats"):
                    reading = True
            else:
                if line.startswith("Ref") == False:
                    parts = line.split(',')
                    ampID = int(parts[0])
                    ampCount = int(parts[1])
                    templateCount = int(parts[2])
                    discardCount = int(parts[3])
                    
                    referenceAmpliconStats.append((ampID, ampCount, templateCount, discardCount))
    
    pprint(referenceAmpliconStats)
    
    with open(filepath) as mutationFile:
        state = 0 # 0 = Mutations 1 = Translocations 2 = Reference Amplicons
        for line in mutationFile:
            if line[:-1] == "Mutations":
                state = 0
                continue
            if line[:-1] == "Translocations":
                state = 1
                continue

            if line[:-1] == "Reference Amplicon Stats":
                state = 2
                break
            if state == 0:
                processMutationLine(line[:-1], personName) # Remove the \n
            elif state == 1:
                processTranslocationLine(line[:-1], personName)
            else:
                processReferenceLine(line[:-1])
                
def processMutationLine(line, personName):
    parts = line.split(", ")
    mutationHash = parts[1]
    ampID = int(mutationHash.strip()[:mutationHash.find(' ')])
    totalOccurrences = referenceAmpliconStats[ampID][1]
    readSplitPt = parts[0].find("/")
    occurrences = int(parts[0][:readSplitPt])
    vaFrequency = round(occurrences / totalOccurrences, 2)
    
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

def processTranslocationLine(line, personName):
    pass
    
def processReferenceLine(line):
    if line.startswith("Count") or line.startswith("ampID"):
        return
    
run()