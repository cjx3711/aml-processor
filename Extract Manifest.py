from csv import *
from itertools import *
from math import pow
from concurrent.futures import *

'''
import difflib
def findAlignCoord(l, r):
	return difflib.SequenceMatcher(None, l, r).find_longest_match(0, len(l), 0, len(r))[:1]
'''

numThreads = 5
complement = {"A" : "T", "C" : "G", "G" : "C", "T" : "A"}
linesPerSeq = {"fastq" : 4}
phredList = [
				"K", "J", "I", "H", "G", "F", "E", 
				"D", "C", "B", "A", "@", "?", ">", 
				"=", "<", ";", ":", "9", "8", "7", 
				"6", "5","4", "3", "2", "1", "0", 
				"/", ".", "-", ",", "+", "*", "(", 
				")", "\'", "&", "%", "$", "#", "\"", "!"
				]
qualityDict = {x : 1 - pow(10, -y/10) for x, y in zip(phredList, range(42, -1, -1))}

def simplifyQuality(p):
	if p >= 0.999:
		return "."
	elif p >= 0.99:
		return ","
	elif p >= 0.9:
		return "?"
	else:
		return "×"

def reverseComplement(bases):
    return ''.join(reversed([complement.get(x, "N") for x in bases]))

def nthAndKthLetter(targetString, targetLetter, n, k):
	lastIndex, nIndex = 0, 0
	for i in range(k):
		if i == n: nIndex = lastIndex
		lastIndex = targetString.index(targetLetter, lastIndex) + 1
	return (nIndex, lastIndex - 1)

def writeToCSV(fileName, targetList):
	with open(fileName, "w+", newline = "") as outFile:
		outWriter = writer(outFile)
		outWriter.writerows(targetList)
		outFile.close()

def grouper(iterable, n):
        args = [iter(iterable)] * n
        return zip(*args)

def IDsAndAlleles():
	with open("Raw/Manifest.csv") as csvFile:
		headerList = next(csvFile)
		csvList = list(reader(csvFile))
		orientedList = [[x[2], reverseComplement(x[3])] if x[1] == "+" else [x[3], reverseComplement(x[2])] for x in csvList]
		combinedList = (y[0][:y[0].rfind(y[1][:15])] + y[1] for y in orientedList)
		tileCounter = [(a, len(tuple(b))) for a, b in groupby([z[0][:-1] for z in csvList])]
		alleleList = []
		for c in tileCounter:
			for d in range(c[1]):
				if d == 0:
					tempAllele = next(combinedList)
				else:
					seqToAdd = next(combinedList)
					tempAllele = tempAllele[:tempAllele.rfind(seqToAdd[:20])] + seqToAdd
			alleleList.append(tempAllele)
		# Decide what to output
		outputList = list(zip([x[0] for x in tileCounter], alleleList))
		writeToCSV("Processed/IDs and Combined Alleles.csv", outputList)

def toJ3X(sequenceData):
	coordIndices = nthAndKthLetter(sequenceData[0], ":", 5, 7)
	sequenceIdentifier = sequenceData[0][coordIndices[0]: coordIndices[1] - 2]
	qualityScores = "".join(simplifyQuality(qualityDict[x]) for x in sequenceData[3][:-1])
	baseSeq = "".join("_" if x == "N" else x for x in sequenceData[1][:-1])
	return sequenceIdentifier + "\n" + baseSeq + "\n" + qualityScores + "\n" + "\n"

def fastqToAny(fastqName, inDir, outDir, conversionFunc):
	linesToSkip = linesPerSeq["fastq"]
	if __name__ ==  "__main__":
		with open(inDir + fastqName) as fastqFile:
			with open(outDir + fastqName[:-6] + ".j3x", "w+", newline = "") as outFile:
				with ProcessPoolExecutor(8) as processManager:
					for x in processManager.map(conversionFunc, grouper(fastqFile, 4), chunksize = 1200):
						outFile.write(x)
					outFile.close()

fastqToAny("AD01_S1_L001_R1_001MINITEST.fastq", "Raw/", "Processed/", toJ3X)