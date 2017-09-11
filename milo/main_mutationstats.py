import glob
from j4xUtils import *
from pprint import pprint
from statistics import median
import json
from genomicsUtils import *
from itertools import groupby, dropwhile
from collections import defaultdict

class MainMutationStats:
                
    def __init__(self, references = 'references/Manifest.csv'):
        self.references = references
        self.statsDir = os.path.join("data", "2-paired")
        self.inDir = os.path.join("data", "3-mutations")
        self.outDir = os.path.join("data", "4-mutationstats")
        self.defaultDictValues = lambda: {
                                            'totalReadsExclDiscards' : 0,
                                            'fileOccurrencePerc': 0,
                                            'fileOccurrences': 0,
                                            'numReads': [],
                                            'numReadsStats': [0,0,0], # Median, min, max
                                            'VAFrequency': [],
                                            'VAFStats': [0,0,0], # Median, min, max
                                            'samples': [],
                                            'coordinates': ''
                                         }
        # self.defaultSigDictValues = lambda: {
        #                                         'fileOccurrences': 0,
        #                                         'mutations': defaultdict(list)
        #                                     }
        self.filenameEnd = "_MUTATIONS.j4x"
        self.mutationHumanDict = defaultdict(self.defaultDictValues)
        self.translocationHumanDict = defaultdict(self.defaultDictValues)
        # self.significantAmpIDDict = defaultdict(self.defaultSigDictValues)

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
    
    def readConfig(self, configFile):
        # Reads config file            
        with open('config.json') as config_file:
            config_data = json.load(config_file)
            if 'j4xstats_minSamples' in config_data:
                self.minSamples = config_data['j4xstats_minSamples']

    def run(self):
        self.readConfig('config.json')
        self.process()
    
    def test(self, configFile, statsDir, inDir, outDir):
        self.statsDir = statsDir
        self.inDir = inDir
        self.outDir = outDir
        self.readConfig(configFile)
        self.process()
        
    def process(self):
        # Reads j3xstats and j4x files
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
        self.processDictDataAnnovar(self.mutationHumanDict, 'annovarStats.csv')
        
        # JJ's code to visualise some stuff
        # sigAmpIDTupleList = [x for x in list(self.significantAmpIDDict.items()) if len([mutation for mutation, samples in x[1]['mutations'].items() if len(samples) < 5]) > 3 ]
        # pprint(sigAmpIDTupleList)

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
            
        outFile.close()
            
    def processSample(self, filepath, statsfilepath, personName):
        self.referenceAmpliconStats = []
        # Calculates discard stats across multiple samples
        with open(statsfilepath) as statsFile:
            refStats = dropwhile(lambda x: x.startswith("Reference Amplicon Stats"), statsFile) # Skip merge and compression statisics
            next(refStats) # Skip info line
            for line in refStats:
                ampID, ampCount, templateCount, discardCount, referenceSequence = line.split(',')
                ampID, ampCount, templateCount, discardCount = (int(num) for num in [ampID, ampCount, templateCount, discardCount]) # Convert strings to int
                
                if len(self.totalReferenceAmpliconStats) <= ampID:
                    self.totalReferenceAmpliconStats.append((ampID, [ampCount], [templateCount], [discardCount], referenceSequence))
                else:
                    self.totalReferenceAmpliconStats[ampID][1].append(ampCount)
                    self.totalReferenceAmpliconStats[ampID][2].append(templateCount)
                    self.totalReferenceAmpliconStats[ampID][3].append(discardCount)
                    
                self.referenceAmpliconStats.append((ampID, ampCount, templateCount, discardCount))

        with open(filepath) as mutationFile:
            # Divide j4x into 3 different parts
            mutations, translocations, refStats = ([line.rstrip() for line in group] for delimiter, group in groupby(mutationFile, key = lambda x: x.rstrip() in ("Mutations", "Translocations", "Reference Amplicon Stats")) if not delimiter)
            for line in mutations:
                    self.processMutationLine(line, personName)
            for line in translocations:
                    self.processTranslocationLine(line, personName)
                    
    def processMutationLine(self, line, personName):
        # Format: "123 / 456", "78 S:91:A-C", "101112" 
        VAFFraction, IDAndMutations, mutationCoords = line.split(", ")
        ampID, mutations = IDAndMutations.split(" ", 1)
        ampID = int(ampID)
        totalReadsInclDiscards = self.referenceAmpliconStats[ampID][1]
        numReads, totalReadsExclDiscards = [int(x.strip()) for x in VAFFraction.split("/")]
        VAFrequency = round(numReads / totalReadsInclDiscards, 2)
        
        self.mutationHumanDict[IDAndMutations]['totalReadsExclDiscards'] += numReads
        self.mutationHumanDict[IDAndMutations]['fileOccurrences'] += 1
        self.mutationHumanDict[IDAndMutations]['numReads'].append(numReads)
        self.mutationHumanDict[IDAndMutations]['VAFrequency'].append(VAFrequency)
        self.mutationHumanDict[IDAndMutations]['samples'].append(personName)
        if ( self.mutationHumanDict[IDAndMutations]['coordinates'] == '' ):
            self.mutationHumanDict[IDAndMutations]['coordinates'] = mutationCoords

        # if self.mutationIsLarge(mutations) or self.mutationIsComplex(mutations):
        #     if VAFrequency > 0.03:
        #         self.significantAmpIDDict[ampID]['fileOccurrences'] += 1
        #         self.significantAmpIDDict[ampID]['mutations'][mutations].append((personName, VAFrequency, numReads))

    # def mutationIsLarge(self, mutations):
    #     return any([len(mutation.split(":")[-1]) > 5 for mutation in mutations.split(" ")])
    # 
    # def mutationIsComplex(self, mutations):
    #     return len(mutations.split(" ")) > 2

    def processTranslocationLine(self, line, personName):
        parts = line.split(", ")
        translocationDescriptor = ', '.join(parts[1:])
        ampID = int(translocationDescriptor.strip()[:translocationDescriptor.find(' ')])
        totalReadsInclDiscards = self.referenceAmpliconStats[ampID][1]
        numReads = int(parts[0])
        vaFrequency = round(numReads / totalReadsInclDiscards, 2)
        
        self.translocationHumanDict[translocationDescriptor]['totalReadsExclDiscards'] += numReads
        self.translocationHumanDict[translocationDescriptor]['fileOccurrences'] += 1
        self.translocationHumanDict[translocationDescriptor]['numReads'].append(numReads)
        self.translocationHumanDict[translocationDescriptor]['VAFrequency'].append(vaFrequency)
        self.translocationHumanDict[translocationDescriptor]['samples'].append(personName)
        
    def createOutFile(self, outputFile):
        if not os.path.exists(self.outDir):
            os.makedirs(self.outDir)
        return open(os.path.join(self.outDir, outputFile), "w+", newline = "")
        
    def processDictDataAnnovar(self, humanDict, outputFile):
        significantTupleList = self.filterList(humanDict)
        outFile = self.createOutFile(outputFile)
        
        for mutID, x in enumerate(significantTupleList):
            mutationHash = x[0]
            coordinates = x[1]['coordinates']
            sampleList = x[1]['samples']
            
            ampID = int(mutationHash.split(' ')[0])
            coordinates = coordinates.split(' ')
            mutations = hashToMutationArray(mutationHash)
            
            if len(mutations) != len(coordinates):
                print('Error: Mutations and Coordinate parts not same length')
            
            for mutation, coordinate in zip(mutations, coordinates):
                chromosome = self.ampliconRefs[ampID][1]
                
                startCoord = int(coordinate)
                endCoord = int(coordinate)
                    
                original = mutation['from'] if len(mutation['from']) > 0 else '-'
                mutated = mutation['to'] if len(mutation['to']) > 0 else '-'
                
                files = 'files:, {0}, {1}%'.format(x[1]['fileOccurrences'], x[1]['fileOccurrencePerc'])
                numReads = 'numReads, {0}, {1}, {2}'.format(x[1]['numReadsStats'][0], x[1]['numReadsStats'][1], x[1]['numReadsStats'][2])
                vaf = 'VAF, {0}, {1}, {2}'.format(x[1]['VAFStats'][0], x[1]['VAFStats'][1], x[1]['VAFStats'][2])
                sampleString = 'Samples, {0}'.format(';'.join(sampleList))
                comments = "comments:, MID:, {0}, AID:, {1}, {2}, {3}, {4}, {5}, {6}".format(mutID, ampID, self.ampliconRefs[ampID][0], files, numReads, vaf, sampleString)
                                
                if mutation['type'] == 'S' and len(original) != len(mutated): # If the length is not the same, split into two different mutations
                    # Addition
                    pwrite(outFile,'{0}   {1}   {2}   {3}   {4}   {5}'.format(chromosome, startCoord, endCoord, '-', mutated, comments ), False)
                    # Deletion
                    endCoord += len(mutation['from']) - 1
                    pwrite(outFile,'{0}   {1}   {2}   {3}   {4}   {5}'.format(chromosome, startCoord, endCoord, original, '-', comments ), False)
                else:
                    if ( mutation['type'] == 'D' ):
                        endCoord += len(mutation['from']) - 1
                    pwrite(outFile,'{0}   {1}   {2}   {3}   {4}   {5}'.format(chromosome, startCoord, endCoord, original, mutated, comments ), False)
        outFile.close()
        
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
            
        outFile.close()

    def filterList(self, humanDict):

        # Format: (ampID + mutation, {VAFrequency, VAFStats, fileOccurrencePerc, fileOccurrences, samples, numReads, numReadsStats, totalReadsInclDiscards})
        mutationTupleList = list(humanDict.items())
        filteredTupleList = [x for x in mutationTupleList if x[1]['fileOccurrences'] >= self.minSamples]
        
        # Stably sort the list by file occurrences,
        filteredTupleList.sort(key = lambda tup: tup[1]['fileOccurrences'], reverse = True)
        # Average VAFrequency
        filteredTupleList.sort(key = lambda tup: sum(tup[1]['VAFrequency']) / tup[1]['fileOccurrences'], reverse = True)
        # 
        filteredTupleList.sort(key = lambda tup: tup[0][:tup[0].find(" ")])
        
        """# If the mutant has a high VAF, or has the same amplicon as a mutant with high VAF and passes a lower VAF threshold
                                significantTupleList = []
                                for x in filteredTupleList:
                                    if (
                                        sum(x[1]['VAFrequency']) / x[1]['fileOccurrences'] > self.VAFThreshold or
                                            (
                                            len(significantTupleList) > 1 and x[0][:x[0].find(" ")] == significantTupleList[-1][0][:significantTupleList[-1][0].find(" ")] and
                                            sum(x[1]['VAFrequency']) / x[1]['fileOccurrences'] > self.VAFThresholdSubordinate
                                            )
                                        ):
                                        significantTupleList.append(x)"""
        return filteredTupleList
        
if __name__ ==  "__main__":
    mutationStats = MainMutationStats()
    mutationStats.run()