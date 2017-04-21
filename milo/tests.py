import unittest

from genomicsUtils import *


# Here's our "unit tests".
class IsOddTests(unittest.TestCase):

    def testOne(self):
        self.assertTrue(IsOdd(1))

    def testTwo(self):
        self.assertTrue(IsOdd(2))

def main():
    unittest.main()

if __name__ == '__main__':
    main()
