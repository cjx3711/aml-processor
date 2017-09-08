"""
This file takes the annovar input file and the annovar output file and combines them into one CSV.
annovarInput is the output file from the j4x stats.
annovarOutput is the output file from the annovar program itself.
"""
import sys
import os


# Step 1: Copy the output from mutationstats to annovar folder


# Step 2: Run the perl annovar function


# Step 3: Take both file outputs and merge them.

outDir = os.path.join("data", "4-mutationstats")

j4xMutationList = []
discoveredMutationList = []
# Read files from program parameters
if len(sys.argv) != 3:
    print('Please provide two arguments.')
else:
    annovarInput = sys.argv[1]
    annovarOutput = sys.argv[2]
    
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