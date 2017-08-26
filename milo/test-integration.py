import unittest
import tempfile
import sys
import os


def compare_files(self, outfile_path, testFile_path):
    with open(outfile_path) as outFile, open(testFile_path) as testFile:
        zipped = zip(outFile.readlines(), testFile.readlines())
        
        for pair in enumerate(zipped):
            try:
                self.assertEqual(pair[1][0], pair[1][1])
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
        pass
    

def main():
    unittest.main()

if __name__ == '__main__':
    main()