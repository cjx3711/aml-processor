import pygtrie
import time

kgramLength = 10
spacing = 10

# indelTree = pygtrie.CharTrie()
indelHash = {}

start_time = time.time()

referenceCount = 0;
with open("references/Manifest.csv", "r", newline = "") as references:

    # Read through each line of the reference file and create the kgrams
    lineno = -1
    for line in references:
        if lineno >= 0:
            referenceCount += 1
            sequence = line.split(',')[2]
            for i in range(0,len(sequence) - kgramLength, spacing):
                kgram = sequence[i:i+kgramLength]
                # # Trie
                # if indelTree.has_key(kgram) == False:
                #     indelTree[kgram] = []
                # indelTree[kgram].append((lineno,i))
                if (kgram in indelHash) == False:
                    indelHash[kgram] = []
                indelHash[kgram].append((lineno,i))
        lineno += 1

    # Add each of the kgrams to the trie


# start_time = time.time()

# with open("data/MINITEST_AD01_S1_L001PAIRED.j3x") as test:
#     linecounter = -1
#     for line in test:
#         linecounter += 1
#         if linecounter % 4 != 2:
#             pass
#
#         matches = []
#         for i in range(0,len(line) - kgramLength + 1):
#             kgram = line[i:i+kgramLength]
#             if indelTree.has_key(kgram):
#                 matches.append(indelTree[kgram])


# print("--- %s seconds ---" % (time.time() - start_time))


badCount = 0
total = 0

with open("data/MINITEST_AD01_S1_L001PAIRED.j3x") as test:
    linecounter = -1
    for line in test:
        linecounter += 1
        if linecounter % 4 != 1:
            continue

        total += 1
        # if linecounter > 2000:
        #     break

        matches = []
        matchCounts = [[0,x] for x in range(referenceCount)]
        
        i = 0
        while i < len(line) - kgramLength + 1:
        # for i in range(0,len(line) - kgramLength + 1):
            kgram = line[i:i+kgramLength]
            if kgram in indelHash:
                currentMatches = indelHash[kgram]
                i += kgramLength
                for match in currentMatches:
                    matchCounts[match[0]][0] += 1
                matches.extend(currentMatches)
            i += 1
        matchCounts.sort(reverse = True)
        if ( matchCounts[1][0] / matchCounts[0][0] > 0.5 ):
            badCount += 1
            # print(str(matchCounts[0]) + str(matchCounts[1]))
        #

        # print(matchCounts)
        # print(matches)

print(badCount)
print(total)

print("--- %s seconds ---" % (time.time() - start_time))
