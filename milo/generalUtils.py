"""
Contains other utils for the application
"""

from DataTypes import *
import os
import sys
import json

def getFileList(fileListJson):
    with open(fileListJson) as configFile:
        filesConfig = json.load(configFile)
        if ( 'process' in filesConfig ):
            process = filesConfig['process']
            return [ generateFilenames(x) for x in process ]
        else:
            print ('Error: files config is in incorrect format')
            return []
    
def generateFilenames(baseName):
    baseName = baseName.strip('_')
    fastq1 = baseName + '_R1_001.fastq'
    fastq2 = baseName + '_R2_001.fastq'
    paired = baseName + '_PAIRED.j3x'
    pairedStats = baseName + '_PAIRED.j3x.stats'
    tiled = baseName + '_TILED.j3x'
    mutation = baseName + '_MUTATIONS.j4x'
    
    return FileTypes(fastq1, fastq2, paired, pairedStats, tiled, mutation)


def pwrite(file, message, shouldPrint = True):
    if shouldPrint:
        print(message)
    file.write(message + "\n")
    
def removeFile(file):
    try:
        os.remove(file)
    except OSError:
        pass

out = None
# Disable printing
def blockPrint():
    global out
    sys.stdout = out = open(os.devnull, 'w')

# Restore
def enablePrint():
    global out
    if ( out != None ): out.close()
    sys.stdout = sys.__stdout__