# from Bio.pairwise2 import *
from j3xUtils import *
from concurrent.futures import *

numThreads = 12

def alignUnpaired(left, right):
    print(left)
    print(right)

    # Slide the logical way
    highestScore = -99999
    finalOverlap = 0
    for overlap in range(5, len(left)):
        finalOverlap = overlap
        score = 0
        for i in zip(left[-overlap:], right[:overlap]):
            if i[0] == 'N' or i[1] == 'N':
                continue
            elif i[0] != i[1]:
                score = score - 5
            elif i[0] == i[1]:
                score = score + 1
            else:
                print("Hello?? " + item)
        # print(score)
        if ( score > highestScore ):
            highestScore = score
        if ( score > 10 ):
            break

    if overlap <= 10:
        for overlap in range(5, len(left)):
            finalOverlap = overlap
            score = 0
            for i in zip(left[:overlap], right[-overlap:]):
                if i[0] == 'N' or i[1] == 'N':
                    continue
                elif i[0] != i[1]:
                    score = score - 5
                elif i[0] == i[1]:
                    score = score + 1
                else:
                    print("Hello?? " + item)
            # print(score)
            if ( score > highestScore ):
                highestScore = score
            if ( score > 10 ):
                break

    print(finalOverlap)
    print(left.ljust(len(left) - finalOverlap + len(left), '-'))
    print(right.rjust(len(right) - finalOverlap + len(right), '-'))
    # return align.localms(left, right, 1, -1, -10, -0.1, one_alignment_only = True)

def pairAligned(aligned):
    baseScore = {"A": 1, "T": 1, "C": 1, "G": 1, "N": 0.5, "-": 0}
    def pickBetter(b1, b2):
        if baseScore[b1] == baseScore[b2] and b1 != b2:
            print("siao liao")
            return b1
        elif baseScore[b1] > baseScore[b2]:
            return b1
        else:
            return b2
    return "".join(pickBetter(*x) for x in zip(aligned[0][:2][0], aligned[0][:2][1]))

def alignAndMerge(sequences):
    return pairAligned(alignUnpaired(sequences[0][1][:-1], reverseComplement(sequences[1][1][:-1]))) + "\n"

def pairToJ3X(fq1, fq2, inDir, outDir):
        with open(inDir + fq1) as fq1File, open(inDir + fq2) as fq2File:
            with open(outDir + fq1[:(fq1.find("_R1_"))] + "PAIRED.j3x", "w+", newline = "") as outFile:
                fq1Iter, fq2Iter = grouper(fq1File, 4), grouper(fq2File, 4)
                with ProcessPoolExecutor(numThreads) as processManager:
                    for x in processManager.map(alignAndMerge, release2Iters(fq1Iter, fq2Iter), chunksize = 250):
                        outFile.write(x)
                    outFile.close()

def release2Iters(iter1, iter2):
    while True:
        i, j = next(iter1), next(iter2)
        if not i or not j: break
        yield i, j

left = [[],"NTGGACTGATATGTGATTTATTCTTTCAACAGCCNNNNNNNGANCCANTGANNAACAAGCTCTCANNNNNNNNCTGAAGATAATGACTCACCNNGGGCCACATTNNAACATTGNAANNTNGCTGGGAGCCTGCACCAAGTCAGGTGGGCTC\n"]
right = [[],"GAGAGTGGAGGATTNAAGNCTGATTGAACNNTTTTCACAACCANNNNNGTNCAGTGAAAATCCTCACTCCNGNTNNGTGAGCCNNNNNGNCNTGGTGCNGGCTCCCAGCANGNNNNNNNNGNNCNNNNGTGGCNNNNNNNNNNNNNNTATN\n"]

# print(alignAndMerge([left, right]))


left = [[],"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAA\n"]
right = [[],"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTT\n"]


left = [[],"AAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"]
right = [[],"TTTTTTTTTTTTTTTTGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n"]

print(alignAndMerge([left, right]))
