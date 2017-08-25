import glob
import os
from j4xUtils import *
from pprint import pprint
from statistics import median
import json
from genomicsUtils import pwrite

statsDir = os.path.join("data", "2-paired")
inDir = os.path.join("data", "3-mutations")
outDir = os.path.join("data", "4-mutationstats")


filenameEnd = "_MUTATIONS.j4x"
mutationHumanDictionary = {}
translocationHumanDictionary = {}

# To find average reference amplicon stats across files
totalReferenceAmpliconStats = []

# Stats from the stats file of the reference amplicon
# Format:
# [
# (ampID, Amplicon Count, Template Count, Discard Count, Discard Rate),
# ... x 571
# ]
referenceAmpliconStats = [] 

totalFiles = 0
minFileOccurences = 1
significanceThreshold = 0.2

def run():
    global totalFiles
    global minFileOccurences
    global significanceThreshold
    global totalReferenceAmpliconStats
    
    with open('config.json') as config_file:
        config_data = json.load(config_file)
        if ( 'j4xstats_minFileOccurences' in config_data ):
            minFileOccurences = config_data['j4xstats_minFileOccurences']
        if ( 'j4xstats_significanceThreshold' in config_data ):
            significanceThreshold = config_data['j4xstats_significanceThreshold']

    for filepath in glob.glob(os.path.join(inDir, '*.j4x')):
        # Read all files in input directory
        filenameParts = filepath.split(os.sep)
        filenameOnly = filenameParts[len(filenameParts)-1]
        if (filenameOnly[-len(filenameEnd):] == filenameEnd):
            personName = filenameOnly[:-len(filenameEnd)]
            statsFile = statsDir + os.sep + personName + '_PAIRED.j3x.stats' 
            processSingleHuman(filepath, statsFile, personName)
            totalFiles += 1
    
    processDictData(mutationHumanDictionary, 'mutationStats.txt')
    processDictData(translocationHumanDictionary, 'translocationStats.txt')

    outFile = open(os.path.join(outDir, 'referenceStats.txt'), "w+", newline = "")
    outFile.write('ampID, Amplicon Count, Discard Count, Discard Rate\n')
    for refStats in totalReferenceAmpliconStats:
        ampCount = median(refStats[1])
        discardCount = median(refStats[3])
        outFile.write(str(refStats[0]))
        outFile.write(', ')
        outFile.write(str(ampCount))
        outFile.write(', ')
        outFile.write(str(discardCount))
        outFile.write(', ')
        if ( ampCount == 0 ):
            outFile.write('-')
        else:
            outFile.write(str(int(discardCount * 1000 / ampCount) / 10))
        outFile.write('%\n')
        
        median(refStats[2])
        median(refStats[3])
        
def processDictData(humanDictionary, outputFile):
    global totalFiles
    global minFileOccurences
    global significanceThreshold
        
    mutationTupleList = list(humanDictionary.items())
    filteredTupleList = [x for x in mutationTupleList if x[1]['fileOccurrences'] >= minFileOccurences]
    
    # Stably sort the list by file occurrences,
    filteredTupleList.sort(key = lambda tup: tup[1]['fileOccurrences'], reverse = True)
    # Average VAFrequency
    filteredTupleList.sort(key = lambda tup: sum(tup[1]['VAFrequency']) / tup[1]['fileOccurrences'], reverse = True)
    # 
    filteredTupleList.sort(key = lambda tup: tup[0][:tup[0].find(" ")])
    significantTupleList = []
    for x in filteredTupleList:
        # High VAF
        if ( sum(x[1]['VAFrequency']) / x[1]['fileOccurrences'] > significanceThreshold ):
            significantTupleList.append(x)
        # Same amplicon number as one with high VAF
        elif ( len(significantTupleList) > 1 and x[0][:x[0].find(" ")] == significantTupleList[-1][0][:significantTupleList[-1][0].find(" ")]):
            significantTupleList.append(x)
            
    # pprint(significantTupleList)
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    outFile = open(os.path.join(outDir, outputFile), "w+", newline = "")
    
    
    for x in significantTupleList:
        # Calculate the stats of the stats
        x[1]['fileOccurencePerc'] = int(x[1]['fileOccurrences'] * 1000 / totalFiles) / 10
        x[1]['occurrenceStats'][0] = median(x[1]['occurrences'])
        x[1]['occurrenceStats'][1] = min(x[1]['occurrences'])
        x[1]['occurrenceStats'][2] = max(x[1]['occurrences'])
        
        x[1]['VAFrequencyStats'][0] = median(x[1]['VAFrequency'])
        x[1]['VAFrequencyStats'][1] = min(x[1]['VAFrequency'])
        x[1]['VAFrequencyStats'][2] = max(x[1]['VAFrequency'])
        
        pwrite(outFile,'{0}\nfiles: {1} ({2}%)'.format(x[0], x[1]['fileOccurrences'], x[1]['fileOccurencePerc']), False)
        pwrite(outFile,'Reads (med, min, max): {0} {1} {2}'.format(x[1]['occurrenceStats'][0], x[1]['occurrenceStats'][1], x[1]['occurrenceStats'][2]), False)
        pwrite(outFile,'VAF   (med, min, max): {0} {1} {2}'.format(x[1]['VAFrequencyStats'][0], x[1]['VAFrequencyStats'][1], x[1]['VAFrequencyStats'][2]), False)
        pwrite(outFile, '', False)
        
def processSingleHuman(filepath, statsfilepath, personName):
    global referenceAmpliconStats
    global totalReferenceAmpliconStats
    referenceAmpliconStats = [] 
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
                    referenceSequence = parts[4]
                    
                    if (len(totalReferenceAmpliconStats) <= ampID ):
                        totalReferenceAmpliconStats.append((ampID, [ampCount], [templateCount], [discardCount], referenceSequence))
                    else:
                        totalReferenceAmpliconStats[ampID][1].append(ampCount)
                        totalReferenceAmpliconStats[ampID][2].append(templateCount)
                        totalReferenceAmpliconStats[ampID][3].append(discardCount)
                        
                    referenceAmpliconStats.append((ampID, ampCount, templateCount, discardCount))

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
            'fileOccurencePerc': 0,
            'fileOccurrences': 0,
            'occurrences': [],
            'occurrenceStats': [0,0,0], # Median, min, max
            'VAFrequency': [],
            'VAFrequencyStats': [0,0,0], # Median, min, max
            'humans': []
        }
    
    mutationHumanDictionary[mutationHash]['totalOccurrences'] += occurrences
    mutationHumanDictionary[mutationHash]['fileOccurrences'] += 1
    mutationHumanDictionary[mutationHash]['occurrences'].append(occurrences)
    mutationHumanDictionary[mutationHash]['VAFrequency'].append(vaFrequency)
    mutationHumanDictionary[mutationHash]['humans'].append(personName)

def processTranslocationLine(line, personName):
    parts = line.split(", ")
    translocationDescriptor = ', '.join(parts[1:])
    ampID = int(translocationDescriptor.strip()[:translocationDescriptor.find(' ')])
    totalOccurrences = referenceAmpliconStats[ampID][1]
    occurrences = int(parts[0])
    vaFrequency = round(occurrences / totalOccurrences, 2)
    
    if ( translocationDescriptor not in translocationHumanDictionary ):
        translocationHumanDictionary[translocationDescriptor] = {
            'totalOccurrences' : 0,
            'fileOccurencePerc': 0,
            'fileOccurrences': 0,
            'occurrences': [],
            'occurrenceStats': [0,0,0], # Median, min, max
            'VAFrequency': [],
            'VAFrequencyStats': [0,0,0], # Median, min, max
            'humans': []
        }
    
    translocationHumanDictionary[translocationDescriptor]['totalOccurrences'] += occurrences
    translocationHumanDictionary[translocationDescriptor]['fileOccurrences'] += 1
    translocationHumanDictionary[translocationDescriptor]['occurrences'].append(occurrences)
    translocationHumanDictionary[translocationDescriptor]['VAFrequency'].append(vaFrequency)
    translocationHumanDictionary[translocationDescriptor]['humans'].append(personName)
    
def processReferenceLine(line):
    if line.startswith("Count") or line.startswith("ampID"):
        return
    
run()