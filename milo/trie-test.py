import pygtrie
import time

kgramLength = 4
spacing = 4

indelTree = pygtrie.CharTrie()
indelHash = {}


reference = [
"ATCGACCGTCGATGCACCATAG",
"ATCCGTCGACGATGGTTGGATC",
"TGCTAGCTTTAAACTGACGTTA",
"GATACGTACCGTGGATGAGTGA",
]
# Read through each line of the reference file and create the kgrams
lineno = 0
for line in reference:
    # Iterate through string and get kgrams
    for i in range(0,len(line) - kgramLength, spacing):
        kgram = line[i:i+kgramLength]
        # Trie
        if indelTree.has_key(kgram) == False:
            indelTree[kgram] = []
        indelTree[kgram].append((lineno,i))
        if (kgram in indelHash) == False:
            indelHash[kgram] = []
        indelHash[kgram].append((lineno,i))
    lineno += 1

# Add each of the kgrams to the trie


test = [
    "ATCGACCGTCGATGCACCATAG",
    "ATCCGTCGACGATGGTTGGATC",
    "TGCTAGCTTTAAACTGACGTTA",
    "GATACGTACCGTGGATGAGTGA",
]

start_time = time.time()

for line in test:
    matches = []
    for i in range(0,len(line) - kgramLength + 1):
        kgram = line[i:i+kgramLength]
        if indelTree.has_key(kgram):
            matches.append(indelTree[kgram])

print("--- %s seconds ---" % (time.time() - start_time))

start_time = time.time()

for line in test:
    for i in range(0,len(line) - kgramLength + 1):
        kgram = line[i:i+kgramLength]
        if kgram in indelHash:
            matches.append(indelHash[kgram])

print("--- %s seconds ---" % (time.time() - start_time))
