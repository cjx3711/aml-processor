import unittest
import tempfile
import sys
import shutil

from main_ampliconid import *
from main_mutationid import *
from main_mutationstats import *
from genomicsUtils import *
from generalUtils import *
from DataTypes import *

def compare_files_unordered(self, outfile_path, testFile_path, group):
    enablePrint()
    outFile = open(outfile_path)
    testFile = open(testFile_path)
    outIter = grouper(outFile, group)
    testIter = grouper(testFile, group)
    testList = []
    for testBatch in testIter:
        string = ""
        for line in testBatch:
            string += line.strip() + '\n'
        testList.append(string)
    
    for outBatch in outIter:
        string = ""
        for line in testBatch:
            string += line.strip() + '\n'
        try:
            self.assertTrue(string in testList)    
        except AssertionError as e:
            raise
        
    outFile.close()
    testFile.close()
    return True
            
    
def compare_files(self, outfile_path, testFile_path, ignore_lines = []):
    with open(outfile_path) as outFile, open(testFile_path) as testFile:
        zipped = zip(outFile.readlines(), testFile.readlines())
        
        for pair in enumerate(zipped):
            if (pair[0] in ignore_lines): continue
            try:
                self.assertEqual(pair[1][0].strip(), pair[1][1].strip())
            except AssertionError as e:
                enablePrint()
                print('Line: ' + str(pair[0]))
                raise
                
        outFile.close()
        testFile.close()
            
        return True
        
def test_j3x_generic(self, filename):
    blockPrint()
    ampliconID = MainAmpliconID()
    outDir = 'test/j3x/'
    inDir = 'test/j3x/'
    configFile = 'test/config.json'
    referenceFile = 'test/simple-manifest.csv'
    filenameArray = [FileTypes(
        filename + '_R1_001.fastq',
        filename + '_R2_001.fastq',
        filename + '.j3x',
        filename + '.j3x.stats',
        '', ''
    )]
    ampliconID.test(inDir, outDir, configFile, referenceFile, filenameArray)
    
    skip = [2]
    result = compare_files_unordered(self, 'test/j3x/' + filename + '.j3x', 'test/j3x/' + filename + '_EXPECTED.j3x', 4)
    if result: os.remove('test/j3x/' + filename + '.j3x')
    result = compare_files(self, 'test/j3x/' + filename + '.j3x.stats', 'test/j3x/' + filename + '_EXPECTED.j3x.stats', skip)
    if result: os.remove('test/j3x/' + filename + '.j3x.stats')
    shutil.rmtree('test/j3x/discarded/')
        
    enablePrint()

def test_j4x_generic(self, filename):
    # blockPrint()
    mutationID = MainMutationID()
    outDir = 'test/j4x/'
    inDir = 'test/j4x/'
    configFile = 'test/config.json'
    referenceFile = 'test/simple-manifest.csv'
    filenameArray = [FileTypes(
        '', '', '', '',
        filename + '.j3x', filename + '.j4x'
    )]
    mutationID.test(inDir, outDir, configFile, referenceFile, filenameArray)
    
    result = compare_files(self, inDir + filename + '.j4x', inDir + filename + '_EXPECTED.j4x')
    if result: os.remove(inDir + filename + '.j4x')


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
        test_j3x_generic(self, 'SIMPLE')
        
    def test_multi(self):
        test_j3x_generic(self, 'MULTI')
        
    def test_multi(self):
        test_j3x_generic(self, 'SIMPLE_MUT')
        
class j4xTests(unittest.TestCase):
    def test_simple(self):
        # Tests for a simple case. Only ampID 3 is present
        test_j4x_generic(self, 'SIMPLE')
        
    def test_multi(self):
        # Tests for multiple ampIDs. All are forward facing.
        test_j4x_generic(self, 'MULTI')
        
    def test_reverse(self):
        # Tests if the reversed sequence is detected correctly
        test_j4x_generic(self, 'REVERSE')
        
    def test_reverse_mut(self):
        # Tests if the reversed sequence is detected correctly
        # Checks for mutations as well.
        # TODO: Check if coordinates are correct
        test_j4x_generic(self, 'REVERSE_MUT')

class mutationStatsTests(unittest.TestCase):
    def setUp(self):
        self.mutationStats = MainMutationStats('test/simple-manifest.csv')
        statsTestDir = 'test/stats'
        self.mutationStats.test('test/config.json', statsTestDir, statsTestDir, statsTestDir)
    
    def test_exist(self):
        self.assertTrue(os.path.exists('test/stats/annovarStats.csv'))
        self.assertTrue(os.path.exists('test/stats/mutationStats.txt'))
        self.assertTrue(os.path.exists('test/stats/referenceStats.txt'))
        self.assertTrue(os.path.exists('test/stats/translocationStats.txt'))
        
    def tearDown(self):
        # pass
        removeFile('test/stats/annovarStats.csv')
        removeFile('test/stats/mutationStats.txt')
        removeFile('test/stats/referenceStats.txt')
        removeFile('test/stats/translocationStats.txt')

def main():
    global out
    unittest.main()
    if ( out != None ): out.close()

if __name__ == '__main__':
    main()