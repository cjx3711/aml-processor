import unittest
import tempfile
import sys
import os




def write_lamb(outfile_path):
    with open(outfile_path, 'w') as outfile:
        outfile.write("Mary had a little lamb.\nLittle Lamb")
        outfile.close()

class LambTests(unittest.TestCase):
    def test_lamb_output(self):
        outfile_path = tempfile.mkstemp()[1]
        
        try:
            write_lamb(outfile_path)
            with open(outfile_path) as outFile, open('test/lamb.txt') as testFile:
                zipped = zip(outFile.readlines(), testFile.readlines())
                
                for pair in enumerate(zipped):
                    try:
                        self.assertEqual(pair[1][0], pair[1][1])
                    except AssertionError as e:
                        print('Line: ' + str(pair[0]))
                        raise

                outFile.close()
                testFile.close()
        finally:
            # NOTE: To retain the tempfile if the test fails, remove
            # the try-finally clauses
            os.remove(outfile_path)
            
    
                
def main():
    unittest.main()

if __name__ == '__main__':
    main()