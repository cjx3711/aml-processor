# Generates the files.json from the raw files path
import json
from os import listdir
from os.path import isfile, join
onlyFiles = [f for f in listdir('data/1-raw') if isfile(join('data/1-raw', f))]

fileset = set()
for f in onlyFiles:
    fparts = f.split("_")
    filename = "_".join(fparts[:-2])
    if filename not in fileset:
        fileset.add(filename)

filelist = list(fileset)

filesjson = {
    "process" : filelist,
    "skip": []
}
with open('files.gen.json', 'w') as outfile:
    json.dump(filesjson, outfile, indent=2, separators=(',', ': '))
    
print("Done.\nScanned {0} FASTQ pairs.\nOutput to files.gen.json.\nPlease rename to files.json to use.".format(len(filelist)))