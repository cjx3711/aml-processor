"""
This file takes the annovar input file and the annovar output file and combines them into one CSV.
annovarInput is the output file from the j4x stats.
annovarOutput is the output file from the annovar program itself.
"""
import sys


j4xMutationList = []
discoveredMutationList = []
# Read files from program parameters
if len(sys.argv) != 3:
    print('Please provide two arguments.')
else:
    annovarInput = sys.argv[1]
    annovarOutput = sys.argv[2]
    
    # Open annovar input file
    inputFile = open(configFile)
    
    # Open annovar output file
    outputFile = open(configFile)