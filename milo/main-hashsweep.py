import time

from AmpliconMatcherHashSweep import *


start_time = time.time()

ampliconMatcher = AmpliconMatcherHashSweep("references/Manifest.csv")

# with open("data/AD01_S1_L001PAIRED.j3x") as test:
with open("data/MINITEST_AD01_S1_L001PAIRED.j3x") as test:
    linecounter = -1
    for line in test:
        linecounter += 1
        if linecounter % 4 != 1:
            continue

        amplicon = ampliconMatcher.findAmplicon(line)

print(ampliconMatcher.noCount)
print(ampliconMatcher.ummCount)
print(ampliconMatcher.badCount)
print(ampliconMatcher.total)




print("--- %s seconds ---" % (time.time() - start_time))
