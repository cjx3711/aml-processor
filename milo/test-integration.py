import unittest
import tempfile
import sys
import os
import shutil

from main_ampliconid import *

def compare_files(self, outfile_path, testFile_path, ignore_lines = []):
    with open(outfile_path) as outFile, open(testFile_path) as testFile:
        zipped = zip(outFile.readlines(), testFile.readlines())
        
        for pair in enumerate(zipped):
            if (pair[0] in ignore_lines): continue
            try:
                self.assertEqual(pair[1][0].strip(), pair[1][1].strip())
            except AssertionError as e:
                print('Line: ' + str(pair[0]))
                raise

        outFile.close()
        testFile.close()
        
# Used to test if the file testing works.
def write_lamb(outfile_path):
    with open(outfile_path, 'w') as outfile:
        outfile.write("Mary had a little lamb.\nLittle Lamb")
        outfile.close()

class LambTests(unittest.TestCase):
    def test_lamb_output(self):
        outfile_path = tempfile.mkstemp()[1]
        write_lamb(outfile_path)
        compare_files(self, outfile_path, 'test/lamb.txt')
        os.remove(outfile_path)
        
class j3xTests(unittest.TestCase):
    def test_simple(self):
        ampliconID = MainAmpliconID()
        outDir = 'test/j3x/'
        inDir = 'test/j3x/'
        configFile = 'test/config.json'
        referenceFile = 'test/simple-manifest.csv'
        filenameArray = {
            "fastq1": "SIMPLE_R1_001.fastq",
            "fastq2": "SIMPLE_R2_001.fastq",
            "paired": "SIMPLE.j3x"
          },
        ampliconID.test(inDir, outDir, configFile, referenceFile, filenameArray)
        
        skip = [2]
        compare_files(self, 'test/j3x/SIMPLE.j3x', 'test/j3x/SIMPLE_EXPECTED.j3x')
        compare_files(self, 'test/j3x/SIMPLE.j3x.stats', 'test/j3x/SIMPLE_EXPECTED.j3x.stats', skip)
        os.remove('test/j3x/SIMPLE.j3x')
        os.remove('test/j3x/SIMPLE.j3x.stats')
        shutil.rmtree('test/j3x/discarded/')

def main():
    unittest.main()

if __name__ == '__main__':
    main()