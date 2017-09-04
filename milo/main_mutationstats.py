import glob
import os
from j4xUtils import *
from pprint import pprint
from statistics import median
import json
from genomicsUtils import *
from itertools import groupby

class MutationStats:
    def __init__(self, references = 'references/Manifest.csv'):
        self.references = references
        self.statsDir = os.path.join("data", "2-paired")
        self.inDir = os.path.join("data", "3-mutations")
        self.outDir = os.path.join("data", "4-mutationstats")

        self.filenameEnd = "_MUTATIONS.j4x"
        self.mutationHumanDict = {}
        self.translocationHumanDict = {}
        self.significantAmpIDDict = {}

        # To find average reference amplicon stats across files
        self.totalReferenceAmpliconStats = []

        # Stats from the stats file of the reference amplicon
        # Format:
        # [
        # (ampID, Amplicon Count, Template Count, Discard Count, Discard Rate),
        # ... x 571
        # ]
        self.referenceAmpliconStats = [] 

        self.totalFiles = 0
        self.minSamples = 1
        self.VAFThreshold = 0.2
        self.VAFThresholdSubordinate = 0.05
        
        self.ampliconRefs = []
        with open(self.references) as refFile:
            for line in refFile:
                csvCells = line.split(',')
                name = csvCells[1]
                chromosome = extractChromosomeNumber(name)
                shortName = name[:name.find('.')]
                self.ampliconRefs.append((shortName, chromosome))

    def run(self):
        # Reads config file            
        with open('config.json') as config_file:
            config_data = json.load(config_file)
            if 'j4xstats_self.minSamples' in config_data:
                self.minSamples = config_data['j4xstats_self.minSamples']
            if 'j4xstats_self.VAFThreshold' in config_data:
                self.VAFThreshold = config_data['j4xstats_self.VAFThreshold']
            if 'j4xstats_self.VAFThresholdSubordinate' in config_data:
                self.VAFThresholdSubordinate = config_data['j4xstats_self.VAFThresholdSubordinate']

        # Reads j3xstats and j4 files
        for filepath in glob.glob(os.path.join(self.inDir, '*.j4x')):
            filename = filepath.split(os.sep)[-1]
            if filename.endswith(self.filenameEnd): # Confirms if file in folder is the j4x mutations file
                personName = filename[:-len(self.filenameEnd)]
                statsFile = self.statsDir + os.sep + personName + '_PAIRED.j3x.stats' 
                self.processSample(filepath, statsFile, personName)
                self.totalFiles += 1
        
        # After processSample populates our dictionaries, process and write them to file
        self.processDictData(self.mutationHumanDict, 'mutationStats.txt')
        self.processDictData(self.translocationHumanDict, 'translocationStats.txt')
        self.processDictDataAnnovar(self.mutationHumanDict, 'annovarStats.txt')

        # Prints summary statistics for amplicons across samples. e.g. discards per amplicon
        outFile = open(os.path.join(self.outDir, 'referenceStats.txt'), "w+", newline = "")
        outFile.write('ampID, Amplicon Count, Discard Count, Discard Rate\n')
        for refStats in self.totalReferenceAmpliconStats:
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
            
    def processSample(self, filepath, statsfilepath, personName):
        self.referenceAmpliconStats = []
        # Calculates discard stats across multiple samples
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
                        
                        if (len(self.totalReferenceAmpliconStats) <= ampID ):
                            self.totalReferenceAmpliconStats.append((ampID, [ampCount], [templateCount], [discardCount], referenceSequence))
                        else:
                            self.totalReferenceAmpliconStats[ampID][1].append(ampCount)
                            self.totalReferenceAmpliconStats[ampID][2].append(templateCount)
                            self.totalReferenceAmpliconStats[ampID][3].append(discardCount)
                            
                        self.referenceAmpliconStats.append((ampID, ampCount, templateCount, discardCount))

        with open(filepath) as mutationFile:
            # Divide j4x into 3 different parts
            mutations, translocations, refStats = [[line.rstrip() for line in group] for delimiter, group in groupby(mutationFile, key = lambda x: x.rstrip() in ("Mutations", "Translocations", "Reference Amplicon Stats")) if not delimiter]
            for line in mutations:
                    self.processMutationLine(line, personName)
            for line in translocations:
                    self.processTranslocationLine(line, personName)
                    
    def processMutationLine(self, line, personName):
        parts = line.split(", ")
        mutationHash = parts[1]
        coordinates = parts[2]
        ampID = int(mutationHash.strip()[:mutationHash.find(' ')])
        totalReadsAcrossSamples = self.referenceAmpliconStats[ampID][1]
        numReadsplitPt = parts[0].find("/")
        numReads = int(parts[0][:numReadsplitPt])
        vaFrequency = round(numReads / totalReadsAcrossSamples, 2)
        
        if ( mutationHash not in self.mutationHumanDict ):
            self.mutationHumanDict[mutationHash] = {
                'coordinates' : coordinates,
                'totalReadsAcrossSamples' : 0,
                'fileOccurrencePerc': 0,
                'fileOccurrences': 0,
                'numReads': [],
                'numReadsStats': [0,0,0], # Median, min, max
                'VAFrequency': [],
                'VAFStats': [0,0,0], # Median, min, max
                'samples': []
            }
        
        self.mutationHumanDict[mutationHash]['totalReadsAcrossSamples'] += numReads
        self.mutationHumanDict[mutationHash]['fileOccurrences'] += 1
        self.mutationHumanDict[mutationHash]['numReads'].append(numReads)
        self.mutationHumanDict[mutationHash]['VAFrequency'].append(vaFrequency)
        self.mutationHumanDict[mutationHash]['samples'].append(personName)

    def processTranslocationLine(self, line, personName):
        parts = line.split(", ")
        translocationDescriptor = ', '.join(parts[1:])
        ampID = int(translocationDescriptor.strip()[:translocationDescriptor.find(' ')])
        totalReadsAcrossSamples = self.referenceAmpliconStats[ampID][1]
        numReads = int(parts[0])
        vaFrequency = round(numReads / totalReadsAcrossSamples, 2)
        
        if ( translocationDescriptor not in self.translocationHumanDict ):
            self.translocationHumanDict[translocationDescriptor] = {
                'totalReadsAcrossSamples' : 0,
                'fileOccurrencePerc': 0,
                'fileOccurrences': 0,
                'numReads': [],
                'numReadsStats': [0,0,0], # Median, min, max
                'VAFrequency': [],
                'VAFStats': [0,0,0], # Median, min, max
                'samples': []
            }
        
        self.translocationHumanDict[translocationDescriptor]['totalReadsAcrossSamples'] += numReads
        self.translocationHumanDict[translocationDescriptor]['fileOccurrences'] += 1
        self.translocationHumanDict[translocationDescriptor]['numReads'].append(numReads)
        self.translocationHumanDict[translocationDescriptor]['VAFrequency'].append(vaFrequency)
        self.translocationHumanDict[translocationDescriptor]['samples'].append(personName)
    
    def filterList(self, humanDict):
        mutationTupleList = list(humanDict.items())
        # Format: (ampID + mutation, {VAFrequency, VAFStats, fileOccurrencePerc, fileOccurrences, samples, numReads, numReadsStats, totalReadsAcrossSamples})
        filteredTupleList = [x for x in mutationTupleList if x[1]['fileOccurrences'] >= self.minSamples]
        
        # Stably sort the list by file occurrences,
        filteredTupleList.sort(key = lambda tup: tup[1]['fileOccurrences'], reverse = True)
        # Average VAFrequency
        filteredTupleList.sort(key = lambda tup: sum(tup[1]['VAFrequency']) / tup[1]['fileOccurrences'], reverse = True)
        # 
        filteredTupleList.sort(key = lambda tup: tup[0][:tup[0].find(" ")])
        
        # If the mutant has a high VAF, or has the same amplicon as a mutant with high VAF and passes a lower VAF threshold
        significantTupleList = []
        for x in filteredTupleList:
            if (
                sum(x[1]['VAFrequency']) / x[1]['fileOccurrences'] > self.VAFThreshold or
                    (
                    len(significantTupleList) > 1 and x[0][:x[0].find(" ")] == significantTupleList[-1][0][:significantTupleList[-1][0].find(" ")] and
                    sum(x[1]['VAFrequency']) / x[1]['fileOccurrences'] > self.VAFThresholdSubordinate
                    )
                ):
                significantTupleList.append(x)
        return significantTupleList
        
    def createOutFile(self, outputFile):
        if not os.path.exists(self.outDir):
            os.makedirs(self.outDir)
        return open(os.path.join(self.outDir, outputFile), "w+", newline = "")
        
    def processDictDataAnnovar(self, humanDict, outputFile):
        significantTupleList = self.filterList(humanDict)
        outFile = self.createOutFile(outputFile)
        
        for i, x in enumerate(significantTupleList):
            mutation = x[0]
            coordinates = x[1]['coordinates']
            
            mutationParts = mutation.split(' ')
            coordinateParts = coordinates.split(' ')
            ampID = int(mutationParts[0])
            mutations = hashToMutationArray(mutation)
            if len(mutations) == 1 and len(coordinateParts) == 1:
                chromosome = self.ampliconRefs[ampID][1]
                startCoord = coordinateParts[0]
                endCoord = coordinateParts[0]
                original = mutations[0]['from'] if len(mutations[0]['from']) > 0 else '-'
                mutated = mutations[0]['to'] if len(mutations[0]['to']) > 0 else '-'
                comments = "comments: MID: {0} AID: {1} {2}".format(i, ampID, self.ampliconRefs[ampID][0])
                
                pwrite(outFile,'{0}   {1}   {2}   {3}   {4}   {5}'.format(chromosome, startCoord, endCoord, original, mutated, comments ))

    def processDictData(self, humanDict, outputFile):           
        significantTupleList = self.filterList(humanDict)
        outFile = self.createOutFile(outputFile)
        
        for x in significantTupleList:
            # Calculate the stats of the stats
            x[1]['fileOccurrencePerc'] = int(x[1]['fileOccurrences'] * 1000 / self.totalFiles) / 10
            x[1]['numReadsStats'][0] = median(x[1]['numReads'])
            x[1]['numReadsStats'][1] = min(x[1]['numReads'])
            x[1]['numReadsStats'][2] = max(x[1]['numReads'])
            
            x[1]['VAFStats'][0] = median(x[1]['VAFrequency'])
            x[1]['VAFStats'][1] = min(x[1]['VAFrequency'])
            x[1]['VAFStats'][2] = max(x[1]['VAFrequency'])
            
            pwrite(outFile,'{0}\nfiles: {1} ({2}%)'.format(x[0], x[1]['fileOccurrences'], x[1]['fileOccurrencePerc']), False)
            pwrite(outFile,'numReads (med, min, max): {0} {1} {2}'.format(x[1]['numReadsStats'][0], x[1]['numReadsStats'][1], x[1]['numReadsStats'][2]), False)
            pwrite(outFile,'VAF   (med, min, max): {0} {1} {2}'.format(x[1]['VAFStats'][0], x[1]['VAFStats'][1], x[1]['VAFStats'][2]), False)
            pwrite(outFile, '', False)

if __name__ ==  "__main__":
    mutationStats = MutationStats()
    mutationStats.run()