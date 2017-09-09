"""
This file takes the annovar input file and the annovar output file and combines them into one CSV.
annovarInput is the output file from the j4x stats.
annovarOutput is the output file from the annovar program itself.
"""
import sys
import os
from genomicsUtils import *

outDir = os.path.join("data", "4-mutationstats")

# Step 1: Check if the required files exist

annovarInput = os.path.join(outDir,'annovarStats.csv')
annovarOutput = os.path.join(outDir,'myanno.hg19_multianno.csv')

# Step 2: Run the perl annovar function
os.system('perl ../annovar/table_annovar.pl ' + annovarInput + ' ../annovar/humandb/ -buildver hg19 -out myanno -remove -protocol refGene,avsnp147,clinvar_20170130,cosmic70 -operation gx,f,f,f -nastring . -csvout -polish')

os.rename("myanno.hg19_multianno.csv", annovarOutput)
removeFile('myanno.invalid_input')
removeFile('myanno.refGene.invalid_input')

# Step 3: Take both file outputs and merge them.
j4xMutationList = []
discoveredMutationList = []

# Open annovar input file
inputFile = open(annovarInput)

# Open annovar output file
outputFile = open(annovarOutput)

# Open output file
combinedFile = open( os.path.join(outDir, 'annovarOut.csv'), "w+", newline = "")


titleLine = next(outputFile).strip().split(',')
titleLine = ','.join(filter(len, titleLine))
titleLine += ',' + ','.join(['Mutation ID', 'AmpID', 'Name', 'Files', 'Files (%)', 'numReads (med)', 'numReads (min)', 'numReads (max)', 'VAF (med)', 'VAF (min)', 'VAF (max)'])
combinedFile.write(titleLine)
combinedFile.write('\n')
for inputIter, outputIter in zip(inputFile, outputFile):
    
    inputParts = inputIter.strip().split(',')
    
    inputCombined = ','.join([inputParts[2], inputParts[4], inputParts[5], inputParts[7], inputParts[8], inputParts[10], inputParts[11], inputParts[12], inputParts[14], inputParts[15], inputParts[16]])

    combinedFile.write(outputIter.strip() + ', ' + inputCombined + '\n')

