from genomicsUtils import *

def run():
    geneID = 1
    with open('references/Manifest.csv') as references:
        with open('references/GeneManifest.csv', "w+", newline = "") as outFile:
        
            first = True
            geneAmpliconList = []
            for line in references:
                if (first):
                    first = False
                    continue
                csvCells = line.split(',')
                ampliconID = csvCells[0]
                direction = csvCells[2]
                ampliconName = csvCells[1]
                sequence = csvCells[3]
                start = int(csvCells[4])
                end = int(csvCells[5])
                
                if ( direction == '-' ):
                    sequence = reverseComplement(sequence)
                
                tileNumber = int(ampliconName[ampliconName.rfind("_")+1:])
                ampliconNameWithoutTile = ampliconName[:ampliconName.rfind("_tile_")]
                if ( tileNumber == 1 ):
                    geneStr = combineAmpliconsintoGene(geneAmpliconList)
                    outFile.write("{0},{1}\n".format(geneID, geneStr))
                    geneID += 1
                    geneAmpliconList = []
                    
                geneAmpliconList.append((ampliconID, ampliconNameWithoutTile, sequence, start, end))

            geneStr = combineAmpliconsintoGene(geneAmpliconList) # Combine the last amplicon
            outFile.write("{0},{1}\n".format(geneID, geneStr))
            outFile.close()
            
            
def combineAmpliconsintoGene(geneAmpliconList):
    if ( len(geneAmpliconList) == 0 ):
        return
    
    ampIDList = []
    firstCoord = 0
    lastCoord = 0
    startCoord = 0
    endCoord = 0
    currentCombined = ''
    for amplicon in geneAmpliconList:
        sequence = amplicon[2]
        newStartCoord = amplicon[3]
        newEndCoord = amplicon[4]
        if ( startCoord == 0 and endCoord == 0 ): # First amplicon
            firstCoord = newStartCoord
            currentCombined = sequence
        else:
            currentCombined += sequence[lastCoord - newStartCoord + 1:]
        startCoord = newStartCoord
        endCoord = lastCoord = newEndCoord
        ampIDList.append(amplicon[0])
    
    ampIDs = " ".join(ampIDList)
    amplicon = geneAmpliconList[0]
    # print("{0} {1} {2}".format(len(geneAmpliconList), amplicon[0], amplicon[1]))
    return "{0},{1},{2},{3},{4}".format(ampIDs, amplicon[1], currentCombined, firstCoord, lastCoord)
    
    
    
run()