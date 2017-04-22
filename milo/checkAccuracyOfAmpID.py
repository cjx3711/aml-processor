from genomicsUtils import *
from Levenshtein import distance

def checkPercentWrong(j3xDir, manifestDir, j3x):
    failedIDCount = 0
    mismatch = 0
    numLines = 0
    with open(j3xDir + j3x) as j3xFile, open(manifestDir + "Manifest.csv") as csvFile:
        headerList = next(csvFile)
        refSeqs = [x[2] for x in list(reader(csvFile))]
        for read in grouper(j3xFile, 4):
            ampID = int(read[0][3:6])
            if ampID == 0:
                failedIDCount += 1
            elif ampID == 41 or ampID == 570:
                if findCorrect(read, refSeqs, 41, 570) != ampID: mismatch += 1
            elif ampID == 539 or ampID == 569:
                if findCorrect(read, refSeqs, 539, 569) != ampID: mismatch += 1
            elif ampID == 368 or ampID == 571:
                if findCorrect(read, refSeqs, 368, 571) != ampID: mismatch += 1
            elif ampID == 188 or ampID == 197:
                if findCorrect(read, refSeqs, 188, 197) != ampID: mismatch += 1
            elif ampID == 137 or ampID == 453:
                if findCorrect(read, refSeqs, 137, 453) != ampID: mismatch += 1
            numLines += 1
    return mismatch/numLines, failedIDCount/numLines

def findCorrect(read, refSeqs, i, j):
    i -= 1
    j -= 1
    return min((distance(read[1], refSeqs[i]), i),(distance(read[1], refSeqs[j]), j))[1]

print(checkPercentWrong("data/Processed/", "references/", "AD01_S1_L001PAIRED.j3x"))