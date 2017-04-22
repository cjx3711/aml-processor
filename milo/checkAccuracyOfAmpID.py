from genomicsUtils import *
from collections import Counter

def checkPercentWrong(j3xDir, manifestDir, j3x, j3xRef):
    unknownCount = 0
    collidedReads = 0
    numReads = 0
    errorList = []
    with open(j3xDir + j3x) as j3xFile, open(j3xDir + j3xRef) as j3xRefFile:
        refIter = grouper(j3xRefFile, 4)
        for read in grouper(j3xFile, 4):
            # Obtain collision and unknown match statistics
            idSeq, baseSeq, qualitySeq, blank = read
            ampID = idSeq[3:6]
            collisions = idSeq[10:11]
            if ampID == "000":
                unknownCount += 1
            if collisions != "0":
                collidedReads += 1
            numReads += 1

            # Compare against "correct" file to determine if false positives follow any distribution.
            refRead = next(refIter)
            refAmpID = refRead[0][3:6]
            if ampID != refAmpID:
                errorList.append(ampID)
    errorCounter = Counter(errorList)
    print(len(errorList))
    print(errorCounter)
    print("Unknowns: " + str(unknownCount/numReads) + ", Reads with collisions: " + str(collidedReads/numReads))

checkPercentWrong("data/Processed/", "references/", "MINITEST_AD01_S1_L001PAIRED_JX.j3x", "MINITEST_AD01_S1_L001PAIRED_MatchedUnknown.j3x")