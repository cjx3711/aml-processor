import unittest

from genomicsUtils import *
from j4xUtils import *
from j3xUtils import *
from AmpliconMatcherProbabilistic import *
from AmpliconMatcherHashSweep import *
from ReadPairer import *
from TranslocatedBlockMatcher import *
from pprint import pprint

# Here's our "unit tests".
class GenomicsUtils(unittest.TestCase):
    # Tests must start with 'test'
    def test_standard_reverse_compliment(self):
        self.assertEqual(reverseComplement("ATCG"), "CGAT")

    def test_wrong_character_reverse_compliment(self):
        self.assertEqual(reverseComplement("ATCG0s"), "NNCGAT")

    def test_empty_reverse_compliment(self):
        self.assertEqual(reverseComplement(""), "")

    def test_longer_reverse_compliment(self):
        self.assertEqual(reverseComplement("AAAAAAAACCCCCCCCTTTTTTTTGGGGGGGG"), "CCCCCCCCAAAAAAAAGGGGGGGGTTTTTTTT")
    
    def test_simple_difference_empty(self):
        d = simpleDistance("", "")
        self.assertEqual(d, 0)
        
    def test_simple_difference_same(self):
        d = simpleDistance("AAAATTTT", "AAAATTTT")
        self.assertEqual(d, 0)
        
    def test_simple_difference_different(self):
        d = simpleDistance("ATATAT", "TATATA")
        self.assertEqual(d, 2)
        
    def test_simple_difference_sub(self):
        d = simpleDistance("ATATTT", "ATATAT")
        self.assertEqual(d, 1)
        
    def test_simple_difference_add(self):
        d = simpleDistance("CATATAT", "ATATAT")
        self.assertEqual(d, 2)
    
    def test_chromosome_number_1(self):
        n = extractChromosomeNumber('NOTCH1.exon.34.line.210.chr9.139390523.139392010_tile_7')
        self.assertEqual(n, '9')
        
    def test_chromosome_number_2(self):
        n = extractChromosomeNumber('SMC3.p.1001.line.349.chr10.112361782.112361834_tile_1')
        self.assertEqual(n, '10')
        
    def test_chromosome_number_3(self):
        n = extractChromosomeNumber('HRAS.exon.3.line.10.chr11.533765.533944_tile_2')
        self.assertEqual(n, '11')

    def test_chromosome_number_4(self):
        n = extractChromosomeNumber('RUNX1.CDS.2.line.89.chr21.36265221.36265260_tile_1')
        self.assertEqual(n, '21')
        
    def test_chromosome_number_5(self):
        n = extractChromosomeNumber('ZRSR2.CDS.11.line.224.chrX.15840853.15841365_tile_1')
        self.assertEqual(n, 'X')
        
class AmpliconMatcherProbabilisticTests(unittest.TestCase):
    def test_file_read(self):
        # Count the amplicons that were read
        ampMat = AmpliconMatcherProbabilistic("references/Manifest.csv")
        self.assertEqual(len(ampMat.refSeqs), 568)
        self.assertEqual(len(ampMat.ampDict1), 568)
        self.assertEqual(len(ampMat.ampDict2), 1704)
        self.assertEqual(len(ampMat.refSeqsTrunc), 568)

class AmpliconMatcherHashSweepTests(unittest.TestCase):
    def test_file_read(self):
        # Count the amplicons that were read
        ampMat = AmpliconMatcherHashSweep("references/Manifest.csv")
        self.assertEqual(ampMat.referenceCount, 568)

    def test_simple_data(self):
        sample_data = [
            'AAAATTTT',
            'CCCCGGGG'
        ]
        ampMat = AmpliconMatcherHashSweep(sample_data, 4, 4)
        self.assertEqual(ampMat.referenceCount, 2)

        # Check if the hash structure is correct
        self.assertEqual(len(ampMat.ampliconRefs), 4)
        for seq in ['AAAA', 'TTTT', 'GGGG', 'CCCC']:
            self.assertTrue(seq in ampMat.ampliconRefs)
        # (ampID, index) ampID is 1-indexed
        self.assertTrue((1, 0) in ampMat.ampliconRefs['AAAA'])
        self.assertTrue((1, 4) in ampMat.ampliconRefs['TTTT'])
        self.assertTrue((2, 0) in ampMat.ampliconRefs['CCCC'])
        self.assertTrue((2, 4) in ampMat.ampliconRefs['GGGG'])

    def test_overlapping_data(self):
        sample_data = [
            'AAAATTTT',
            'TTTTGGGG'
        ]
        ampMat = AmpliconMatcherHashSweep(sample_data, 4, 4)
        self.assertEqual(ampMat.referenceCount, 2)
        self.assertEqual(len(ampMat.ampliconRefs), 3)
        
        # Check if the hash structure is correct
        for seq in ['AAAA', 'TTTT', 'GGGG']:
            self.assertTrue(seq in ampMat.ampliconRefs)
        self.assertFalse('CCCC' in ampMat.ampliconRefs)
        self.assertEqual(len(ampMat.ampliconRefs['TTTT']), 2)
        # (ampID, index) ampID is 1-indexed
        self.assertTrue((1, 0) in ampMat.ampliconRefs['AAAA'])
        self.assertTrue((1, 4) in ampMat.ampliconRefs['TTTT'])
        self.assertTrue((2, 0) in ampMat.ampliconRefs['TTTT'])
        self.assertTrue((2, 4) in ampMat.ampliconRefs['GGGG'])

    def test_skip_data(self):
        sample_data = [
            'AAAATTTTCCCCGGGG'
        ]
        ampMat = AmpliconMatcherHashSweep(sample_data, 4, 5)
        self.assertEqual(ampMat.referenceCount, 1)

        # Check if the hash structure is correct
        self.assertEqual(len(ampMat.ampliconRefs), 3)
        for seq in ['AAAA', 'TTTC', 'CCGG']:
            self.assertTrue(seq in ampMat.ampliconRefs)
        # (ampID, index) ampID is 1-indexed
        self.assertTrue((1, 0) in ampMat.ampliconRefs['AAAA'])
        self.assertTrue((1, 5) in ampMat.ampliconRefs['TTTC'])
        self.assertTrue((1, 10) in ampMat.ampliconRefs['CCGG'])

    def test_scoring(self):
        sample_data = [
            'AAAATTTTCCCC',
            'CCAATTTTGGGG'
        ]
        amplicon = "AATTTTCCCC"
        ampMat = AmpliconMatcherHashSweep(sample_data, 4, 4)
        self.assertEqual(ampMat.referenceCount, 2)

        # Check if the hash structure is correct
        for seq in ['AAAA', 'TTTT', 'CCCC', 'CCAA', 'GGGG']:
            self.assertTrue(seq in ampMat.ampliconRefs)
            
        largest1, largest2 = ampMat.calculateMatchScore(amplicon)
        self.assertEqual([4, 1], largest1)
        self.assertEqual([2, 2], largest2)
        
class miscj3xUtils(unittest.TestCase):
    def extract_0_id(self):
        ampID1, ampID2 = extractAmpID('ID:000')
        self.assertEqual(ampID1, 0)
        self.assertEqual(ampID2, 0)
        
    def extract_001_id(self):
        ampID1, ampID2 = extractAmpID('ID:001')
        self.assertEqual(ampID1, 1)
        self.assertEqual(ampID2, 0)
        
    def extract_014_id(self):
        ampID1, ampID2 = extractAmpID('ID:014')
        self.assertEqual(ampID1, 14)
        self.assertEqual(ampID2, 0)
    
    def extract_241_id(self):
        ampID1, ampID2 = extractAmpID('ID:241')
        self.assertEqual(ampID1, 241)
        self.assertEqual(ampID2, 0)
        
    def extract_241_111_tl(self):
        ampID1, ampID2 = extractAmpID('TL:241/111')
        self.assertEqual(ampID1, 241)
        self.assertEqual(ampID2, 111)
        
    def extract_121_001_tl(self):
        ampID1, ampID2 = extractAmpID('TL:121/001')
        self.assertEqual(ampID1, 121)
        self.assertEqual(ampID2, 1)
        
class miscj4xUtils(unittest.TestCase):
    def weirdly_specific_function_1(self):
        coords = convertHashPositionsToCoordinates('291 S:29:T-A S:102:T-G', 0)
        self.assertEqual(coords, '29 102')
        
    def weirdly_specific_function_2(self):
        coords = convertHashPositionsToCoordinates('291 S:29:T-A S:102:T-G', 100)
        self.assertEqual(coords, '129 202')
        
class PairedFASTQAligner(unittest.TestCase):
    
    def setUp(self):
        self.readPairer = ReadPairer('test/config.json')
    
    def test_scoring(self):
        left = "AAAAAAAAAAAAA"
        right = "AAAAAAAAAAAAA"
        overlapPairs = tuple(zip(left, right))
        score = self.readPairer.calcScore(overlapPairs, len(left), 0, False)
        self.assertEqual(score, 13)
    
    def test_score_ignore(self):
        left = "AAAAAAAAAAAAA"
        right = "CCCCCCCCCCCCC"
        overlapPairs = tuple(zip(left, right))
        score = self.readPairer.calcScore(overlapPairs, len(left), 0, False)
        self.assertEqual(score, -5)
    
    def test_score_mismatch(self):
        left = "AAAAAAAAAAAAA"
        right = "AAAAAAAAAAAAC"
        overlapPairs = tuple(zip(left, right))
        score = self.readPairer.calcScore(overlapPairs, len(left), 0, False)
        self.assertEqual(score, 7)
    
    def test_no_overlap(self):
        left = "A" * 151
        right = "C" * 151
        lQuality = "K" * len(left)
        rQuality = "K" * len(right)
    
        mergedSequence, qualityScores, noOfCol, newScore = self.readPairer.mergeUnpaired(left, right, lQuality, rQuality, False)
        expected = left + " " + right
    
        # No overlap because too small
        self.assertEqual(mergedSequence, expected)
        self.assertEqual(noOfCol, "?")
    
    def test_small_overlap(self):
        overlap = 5
        sides = 146
        left = "A" * sides + "C" * overlap
        right = "C" * overlap + "G" * sides
        lQuality = "K" * len(left)
        rQuality = "K" * len(right)
    
        mergedSequence, qualityScores, noOfCol, newScore = self.readPairer.mergeUnpaired(left, right, lQuality, rQuality, False)
    
        expected = left + " " + right
    
        # No overlap because too small
        self.assertEqual(mergedSequence, expected)
        self.assertEqual(noOfCol, "?")
    
    def test_normal_overlap(self):
        overlap = 11
        sides = 140
        left = "A" * sides + "C" * overlap
        right = "C" * overlap + "G" * sides
        lQuality = "K" * len(left)
        rQuality = "K" * len(right)
    
        mergedSequence, qualityScores, noOfCol, newScore = self.readPairer.mergeUnpaired(left, right, lQuality, rQuality, False)
    
        expected = "A" * sides + "C" * overlap + "G" * sides
    
        self.assertEqual(mergedSequence, expected)
    
    def test_large_overlap(self):
        overlap = 51
        sides = 100
        left = "A" * sides + "C" * overlap
        right = "C" * overlap + "G" * sides
        lQuality = "K" * len(left)
        rQuality = "K" * len(right)
    
        mergedSequence, qualityScores, noOfCol, newScore = self.readPairer.mergeUnpaired(left, right, lQuality, rQuality, False)
    
        expected = "A" * sides + "C" * overlap + "G" * sides
    
        # vulnerable to homopolymer sequences. It should give the wrong value if the overlap is too big
        self.assertNotEqual(mergedSequence, expected)
    
    def test_actual_sequence(self):
        left = "CAAGTTGATGGGGGCAAAACCAAAATAATTTTCATTTTTAATATACCACACAACACATCTATCTACAAATGCTTTACATTAAATTTACCTGCCAACTGTTTAGCCTGGCTTGCGTTTTCAGTTTGTTTTGTACGTGATGGGGCTGACTTTT"
        right = reverseComplement("AGTCAGGATGTTAGCAGAGCCAGTCAAGACTTGCCGACAAAGGAAACTAGAAGCCAAGAAAGCTGCAGCTGAAAAGCTTTCCTCCCTGGAGAACAGCTCAAATAAAAATGAAAAGGAAAAGTCAGCCCCATCACGTACAAAACAAACTGAA")
        lQuality = "CCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
        rQuality = "CCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGDGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGDFG"[::-1]
        mergedSequence, qualityScores, noOfCol, newScore = self.readPairer.mergeUnpaired(left, right, lQuality, rQuality, False)
    
        expectedMerged = "CAAGTTGATGGGGGCAAAACCAAAATAATTTTCATTTTTAATATACCACACAACACATCTATCTACAAATGCTTTACATTAAATTTACCTGCCAACTGTTTAGCCTGGCTTGCGTTTTCAGTTTGTTTTGTACGTGATGGGGCTGACTTTTCCTTTTCATTTTTATTTGAGCTGTTCTCCAGGGAGGAAAGCTTTTCAGCTGCAGCTTTCTTGGCTTCTAGTTTCCTTTGTCGGCAAGTCTTGACTGGCTCTGCTAACATCCTGACT"
        expecedQuality = "." * len(expectedMerged)
        expectedCol = '0'
        self.assertEqual(mergedSequence, expectedMerged)
        self.assertEqual(qualityScores, expecedQuality)
        self.assertEqual(noOfCol, expectedCol)
    
    def test_actual_sequence_2(self):
        left = "TCTTTCTGCCTCATGCTCTCTCCAACAGGCTTGCAGCCAATTTACTGGAGCAGGGATGACGTAGCCCAGTGGCTCATGTGGGCTGAAAATGAGTTTTCTTTAAGGCCAATTGACAGCAACACATTTGAAATGAATGGCAAAGCTCTCCTGC"
        right = reverseComplement("GGGGTGTTAAAGACCAACCACTAACTAAGAGATTTTCCAAGTTGGGCATATGCCAAGAGTCCAGACTCTCACCTGAATGAGGAGATCGATAGCGAAAGTCCTCTTTGGTCAGCAGCAGGAGAGCTTTGCCATTCATTTCAACTGTGTTGAT")
        lQuality = "CCCCCGGGFCCEF;F9AGFGFDEGGGGGDFAAEFF<FG<CFFFFGFEG8F@C@C7<C:FC7FD8<FDCAC<FEGGC,CFGACFFGGG<F6C<C9<6C@@FF@CCFFGEGGGGGGGGFG8=?4EEFFFCCCFGGFFF<EGGGGCF8?==BC@"
        rQuality = "-8BC8B@FC<9E--6FFECEFGGG8CCFCGGG,CFFF,,,;EFFCCC,CFF6CFCA8,CC,C96C@EF9,C,,66,<,E<,@FB8DGGEFE@<@C@:>?,9,<?FFGCGGG?,<<E<,9F:<?@<AEG9A?CEF,AFGECF,@EEFEDF,A"[::-1]
    
        mergedSequence, qualityScores, noOfCol, newScore = self.readPairer.mergeUnpaired(left, right, lQuality, rQuality, False)
    
        expectedMerged = "TCTTTCTGCCTCATGCTCTCTCCAACAGGCTTGCAGCCAATTTACTGGAGCAGGGATGACGTAGCCCAGTGGCTCATGTGGGCTGAAAATGAGTTTTCTTTAAGGCCAATTGACAGCAACACATTTGAAATGAATGGCAAAGCTCTCCTGCTGCTGACCAAAGAGGACTTTCGCTATCGATCTCCTCATTCAGGTGAGAGTCTGGACTCTTGGCATATGCCCAACTTGGAAAATCTCTTAGTTAGTGGTTGGTCTTTAACACCCC"
        expecedQuality = ".............,.,...................,..,.........,.....,,.,..,..,,.....,.....?..........,.,.,.,,,......................,,.?..............,.......,.,,...,?.........,?,?.,,...,.......,...?,.?,?,,??.?,....,,.?..?,....,...?......,???....?.......,.........,??.,,....,..,?"
        expectedCol = '0'
        self.assertEqual(mergedSequence, expectedMerged)
        self.assertEqual(qualityScores, expecedQuality)
        self.assertEqual(noOfCol, expectedCol)
    
    def test_actual_sequence_delete(self):
        left = "TGCAGTCCCAGCCCACAGCCCCCCTCCTCCCTCAGACTCAGGAGTCCATAGCGAATTTCGACGATCGTTGCATTAACTCGCGAACAAACGGATCTCGTATGCCGTCTTCTGCTTGAAAACAAAAACAATTTCAATACGTCTGAAAATGTTA"
        right = reverseComplement("TATGGACTCCTGAGTCTGAGGGAGGAGGGGGGCTGTGGGCTGGGACTGCAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAAGGTTCAGTGTAGATCTCGGTGGGCGCCGGATCATTAATAAAAAAAAAACAGAAATTCACGGTGAGC")
        lQuality = "CCCCCGEEFGGGGGGGGEGCCEFGCEGGFCFGGGGDG@EFGFGGGAFGGGEGGGGCDEEGFGGEEEFEG@CFEDEEGCFGCEFE@FEGFGG=4+,EFF7D=FGG=CEFB??,CFGFGG<,?,,4B+:,,,8<5,,C,B:,,+5,:,:,@B,"
        rQuality = "CCCCCFGGGGGGFGGFGFFGGGGGG;CEDEGGGGGGFGGGGGDG<AEFG?EBFFGGG:@CEGCF>FECFFFGGE<BFF,BCE,EABFB:>>FDFF9FF,=3@F<>:C*3<*5<****8,<;<,@,8;:=E*1*1*++2,27,,,***4*2+"[::-1]
    
        mergedSequence, qualityScores, noOfCol, newScore = self.readPairer.mergeUnpaired(left, right, lQuality, rQuality, False)
    
        expectedMerged = "TGCAGTCCCAGCCCACAGCCCCCCTCCTCCCTCAGACTCAGGAGTCCATA"
        expecedQuality = ".................................................."
        expectedCol = '0'
        self.assertEqual(mergedSequence, expectedMerged)
        self.assertEqual(qualityScores, expecedQuality)
        self.assertEqual(noOfCol, expectedCol)

class TranslocatedBlockMatcherTests(unittest.TestCase):
    def setUp(self):
        self.translocatedBlockMatcher = TranslocatedBlockMatcher()

    def test_1(self): # Refs completely matched to read
        read = "ABCDEFGH1234567"
        ref1 = "1234567"
        ref2 = "ABCDEFGH"
        
        matchingBlocks = self.translocatedBlockMatcher.findTranslocatedMatchingBlocks(read, ref1, ref2)
        self.assertEqual(len(matchingBlocks), 2)
        blocks = []
        for match in matchingBlocks:
            blocks.append(str(match))
        
        self.assertTrue("('R2', Match(a=0, b=0, size=8))" in blocks)
        self.assertTrue("('R1', Match(a=8, b=0, size=7))" in blocks)

    def test_2(self): # Refs incompletely matched to read
        read = "ABCDEFGH1234567"
        ref1 = "12345670000000"
        ref2 = "ZZZZZZABCDEFGH"
        
        matchingBlocks = self.translocatedBlockMatcher.findTranslocatedMatchingBlocks(read, ref1, ref2)
        self.assertEqual(len(matchingBlocks), 2)
        blocks = []
        for match in matchingBlocks:
            blocks.append(str(match))
        
        self.assertTrue("('R2', Match(a=0, b=6, size=8))" in blocks)
        self.assertTrue("('R1', Match(a=8, b=0, size=7))" in blocks)

    def test_3(self): # Same as test_2, but with unpaired areas on opposite ends
        read = "ABCDEFGH1234567"
        ref1 = "00001234567"
        ref2 = "ABCDEFGHZZZZZ"
        
        matchingBlocks = self.translocatedBlockMatcher.findTranslocatedMatchingBlocks(read, ref1, ref2)
        self.assertEqual(len(matchingBlocks), 2)
        blocks = []
        for match in matchingBlocks:
            blocks.append(str(match))
        
        self.assertTrue("('R2', Match(a=0, b=0, size=8))" in blocks)
        self.assertTrue("('R1', Match(a=8, b=4, size=7))" in blocks)

    def test_4(self): # 1bp insertions on both refs
        read = "ABCDEFGH1234567"
        ref1 = "12345687"
        ref2 = "ABCZDEFGH"
        
        matchingBlocks = self.translocatedBlockMatcher.findTranslocatedMatchingBlocks(read, ref1, ref2)
        #self.assertEqual(len(matchingBlocks), 4)
        blocks = []
        for match in matchingBlocks:
            blocks.append(str(match))
        
        self.assertTrue("('R1', Match(a=8, b=0, size=6))" in blocks)
        self.assertTrue("('R2', Match(a=3, b=4, size=5))" in blocks)
        self.assertTrue("('R2', Match(a=0, b=0, size=3))" in blocks)
        self.assertTrue("('R1', Match(a=14, b=7, size=1))" in blocks)

    def test_5(self): # Insertions on both refs, with one >1bp long
        read = "ABCDEFGH1234567"
        ref1 = "12345688887"
        ref2 = "ABCZDEFGH"
        
        matchingBlocks = self.translocatedBlockMatcher.findTranslocatedMatchingBlocks(read, ref1, ref2)
        #self.assertEqual(len(matchingBlocks), 4)
        blocks = []
        for match in matchingBlocks:
            blocks.append(str(match))
        
        self.assertTrue("('R1', Match(a=8, b=0, size=6))" in blocks)
        self.assertTrue("('R2', Match(a=3, b=4, size=5))" in blocks)
        self.assertTrue("('R2', Match(a=0, b=0, size=3))" in blocks)
        self.assertTrue("('R1', Match(a=14, b=10, size=1))" in blocks)

    def test_6(self): # 1bp insertiion on one ref
        read = "ABCDEFGH1234567"
        ref1 = "1234567"
        ref2 = "ABCZDEFGH"
        
        matchingBlocks = self.translocatedBlockMatcher.findTranslocatedMatchingBlocks(read, ref1, ref2)
        self.assertEqual(len(matchingBlocks), 3)
        blocks = []
        for match in matchingBlocks:
            blocks.append(str(match))
        
        self.assertTrue("('R1', Match(a=8, b=0, size=7))" in blocks)
        self.assertTrue("('R2', Match(a=3, b=4, size=5))" in blocks)
        self.assertTrue("('R2', Match(a=0, b=0, size=3))" in blocks)

    def test_7(self): # 4bp + 1bp insertion on one ref
        read = "ABCDEFGH1234567"
        ref1 = "1234567"
        ref2 = "ABCZZZZDEFGZH"
        
        matchingBlocks = self.translocatedBlockMatcher.findTranslocatedMatchingBlocks(read, ref1, ref2)
        self.assertEqual(len(matchingBlocks), 4)
        blocks = []
        for match in matchingBlocks:
            blocks.append(str(match))
        
        self.assertTrue("('R1', Match(a=8, b=0, size=7))" in blocks)
        self.assertTrue("('R2', Match(a=3, b=7, size=4))" in blocks)
        self.assertTrue("('R2', Match(a=0, b=0, size=3))" in blocks)
        self.assertTrue("('R2', Match(a=7, b=12, size=1))" in blocks)

    def test_8(self): # Multiple insertions of varying lengths on both refs
        read = "ABCDEFGH1234567"
        ref1 = "12000034560007"
        ref2 = "ABCZZDEFGHZZZZZ"
        
        matchingBlocks = self.translocatedBlockMatcher.findTranslocatedMatchingBlocks(read, ref1, ref2)
        #self.assertEqual(len(matchingBlocks), 5)
        blocks = []
        for match in matchingBlocks:
            blocks.append(str(match))
        
        self.assertTrue("('R2', Match(a=3, b=5, size=5))" in blocks)
        self.assertTrue("('R1', Match(a=10, b=6, size=4))" in blocks)
        self.assertTrue("('R2', Match(a=0, b=0, size=3))" in blocks)
        self.assertTrue("('R1', Match(a=8, b=0, size=2))" in blocks)
        self.assertTrue("('R1', Match(a=14, b=13, size=1))" in blocks)

    def test_9(self): # 1bp deletion on one ref
        read = "ABCDEFGH1234567"
        ref1 = "123467"
        ref2 = "ABCDEFGH"
        
        matchingBlocks = self.translocatedBlockMatcher.findTranslocatedMatchingBlocks(read, ref1, ref2)
        #self.assertEqual(len(matchingBlocks), 3)
        blocks = []
        for match in matchingBlocks:
            blocks.append(str(match))
        
        self.assertTrue("('R2', Match(a=0, b=0, size=8))" in blocks)
        self.assertTrue("('R1', Match(a=8, b=0, size=4))" in blocks)
        self.assertTrue("('R1', Match(a=13, b=4, size=2))" in blocks)

    def test_10(self): # >1bp deletion on one ref
        read = "ABCDEFGH1234567"
        ref1 = "1237"
        ref2 = "ABCDEFGH"
        
        matchingBlocks = self.translocatedBlockMatcher.findTranslocatedMatchingBlocks(read, ref1, ref2)
        #self.assertEqual(len(matchingBlocks), 3)
        blocks = []
        for match in matchingBlocks:
            blocks.append(str(match))
        
        self.assertTrue("('R2', Match(a=0, b=0, size=8))" in blocks)
        self.assertTrue("('R1', Match(a=8, b=0, size=3))" in blocks)
        self.assertTrue("('R1', Match(a=14, b=3, size=1))" in blocks)

    def test_11(self): # 1bp deletion on both refs
        read = "ABCDEFGH1234567"
        ref1 = "123467"
        ref2 = "ACDEFGH"
        
        matchingBlocks = self.translocatedBlockMatcher.findTranslocatedMatchingBlocks(read, ref1, ref2)
        #self.assertEqual(len(matchingBlocks), 4)
        blocks = []
        for match in matchingBlocks:
            blocks.append(str(match))
        
        self.assertTrue("('R2', Match(a=2, b=1, size=6))" in blocks)
        self.assertTrue("('R1', Match(a=8, b=0, size=4))" in blocks)
        self.assertTrue("('R1', Match(a=13, b=4, size=2))" in blocks)
        self.assertTrue("('R2', Match(a=0, b=0, size=1))" in blocks)

    def test_12(self): # >1bp deletion on both refs
        read = "ABCDEFGH1234567"
        ref1 = "1237"
        ref2 = "ABCDEH"
        
        matchingBlocks = self.translocatedBlockMatcher.findTranslocatedMatchingBlocks(read, ref1, ref2)
        #self.assertEqual(len(matchingBlocks), 4)
        blocks = []
        for match in matchingBlocks:
            blocks.append(str(match))
        
        self.assertTrue("('R2', Match(a=0, b=0, size=5))" in blocks)
        self.assertTrue("('R1', Match(a=8, b=0, size=3))" in blocks)
        self.assertTrue("('R1', Match(a=14, b=3, size=1))" in blocks)
        self.assertTrue("('R2', Match(a=7, b=5, size=1))" in blocks)

    def test_13(self): # Multiple deletions of varying lengths on both refs
        read = "ABCDEFGH1234567"
        ref1 = "12347"
        ref2 = "ABDE"
        
        matchingBlocks = self.translocatedBlockMatcher.findTranslocatedMatchingBlocks(read, ref1, ref2)
        #self.assertEqual(len(matchingBlocks), 4)
        blocks = []
        for match in matchingBlocks:
            blocks.append(str(match))
        
        self.assertTrue("('R1', Match(a=8, b=0, size=4))" in blocks)
        self.assertTrue("('R2', Match(a=0, b=0, size=2))" in blocks)
        self.assertTrue("('R2', Match(a=3, b=2, size=2))" in blocks)
        self.assertTrue("('R1', Match(a=14, b=4, size=1))" in blocks)

    def test_14(self): # >1bp deletions on ends
        read = "ABCDEFGH1234567"
        ref1 = "1234"
        ref2 = "CDEFGH"
        
        matchingBlocks = self.translocatedBlockMatcher.findTranslocatedMatchingBlocks(read, ref1, ref2)
        self.assertEqual(len(matchingBlocks), 2)
        blocks = []
        for match in matchingBlocks:
            blocks.append(str(match))
        
        self.assertTrue("('R2', Match(a=2, b=0, size=6))" in blocks)
        self.assertTrue("('R1', Match(a=8, b=0, size=4))" in blocks)

    def test_15(self): # Same as test_14, but with deletions flipped to opposite ends
        read = "ABCDEFGH1234567"
        ref1 = "4567"
        ref2 = "ABCDEF"
        
        matchingBlocks = self.translocatedBlockMatcher.findTranslocatedMatchingBlocks(read, ref1, ref2)
        self.assertEqual(len(matchingBlocks), 2)
        blocks = []
        for match in matchingBlocks:
            blocks.append(str(match))
        
        self.assertTrue("('R2', Match(a=0, b=0, size=6))" in blocks)
        self.assertTrue("('R1', Match(a=11, b=0, size=4))" in blocks)
        
    def test_simple(self):
        read = "ABCZDEFGH1234567809"
        ref1 = "0000000000001234567890"
        ref2 = "ZABCDEFGHIJKLMNOPQR"
        
        matchingBlocks = self.translocatedBlockMatcher.findTranslocatedMatchingBlocks(read, ref1, ref2)
        self.assertEqual(len(matchingBlocks), 4)
        blocks = []
        for match in matchingBlocks:
            blocks.append(str(match))
            
        self.assertTrue("('R1', Match(a=9, b=12, size=8))" in blocks)
        self.assertTrue("('R2', Match(a=4, b=4, size=5))" in blocks)
        self.assertTrue("('R2', Match(a=0, b=1, size=3))" in blocks)
        self.assertTrue("('R1', Match(a=17, b=21, size=1))" in blocks)
        
    def test_end_translocation(self):
        read = "ESTEST000012345TESTESTES"
        ref1 = "0000000000001234567890"
        ref2 = "TESTESTESTESTESTESTEST"
        
        matchingBlocks = self.translocatedBlockMatcher.findTranslocatedMatchingBlocks(read, ref1, ref2)

        blocks = []
        for match in matchingBlocks:
            blocks.append(str(match))
        self.assertTrue("('R2', Match(a=15, b=0, size=9))" in blocks)
        self.assertTrue("('R1', Match(a=6, b=8, size=9))" in blocks)

class MutationApplication(unittest.TestCase):
    # Test when mutant in front
    def test_1(self):
        base = "AAATTTCCCGGG"
        mutation = "I:0:-T"
        mutated = mutateSequenceByHash(base,mutation)
        self.assertEqual(mutated, 'TAAATTTCCCGGG')
        
    def test_2(self):
        base = "AAATTTCCCGGG"
        mutation = "S:0:A-T"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'TAATTTCCCGGG')
        
    def test_3(self):
        base = "AAATTTCCCGGG"
        mutation = "D:0:A-"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'AATTTCCCGGG')

    # Test expectation tracking when mutant in front
    # Front mutant is Ins
    def test_4(self):
        base = "AAATTTCCCGGG"
        mutation = "I:0:-T I:9:-TT"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'TAAATTTCCCTTGGG')
        
    def test_5(self):
        base = "AAATTTCCCGGG"
        mutation = "I:0:-T I:9:-TT I:12:-TTT"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'TAAATTTCCCTTGGGTTT')
        
    def test_6(self):
        base = "AAATTTCCCGGG"
        mutation = "I:0:-T I:9:-TT D:11:G-"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'TAAATTTCCCTTGG')
        
    def test_7(self):
        base = "AAATTTCCCGGG"
        mutation = "I:0:-T D:3:TTT- I:9:-TT"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'TAAACCCTTGGG')
        
    def test_8(self):
        base = "AAATTTCCCGGG"
        mutation = "I:0:-T D:3:TTTC- D:10:GG-"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'TAAACCG')
        
    # Test expectation tracking when mutant in front
    # Front mutant is Sub
    def test_9(self):
        base = "AAATTTCCCGGG"
        mutation = "S:0:A-T I:9:-TT"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'TAATTTCCCTTGGG')
        
    def test_10(self):
        base = "AAATTTCCCGGG"
        mutation = "S:0:A-T I:9:-TT I:12:-TTT"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'TAATTTCCCTTGGGTTT')
        
    def test_11(self):
        base = "AAATTTCCCGGG"
        mutation = "S:0:A-T I:9:-TT D:11:G-"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'TAATTTCCCTTGG')
        
    def test_12(self):
        base = "AAATTTCCCGGG"
        mutation = "I:0:-T D:2:ATTT- I:9:-TT"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'TAACCCTTGGG')
        
    def test_13(self):
        base = "AAATTTCCCGGG"
        mutation = "I:0:-T D:2:ATTTC- D:10:GG-"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'TAACCG')

    # Front mutant is Del
    def test_14(self):
        base = "AAATTTCCCGGG"
        mutation = "D:0:A- I:9:-TT"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'AATTTCCCTTGGG')
        
    def test_15(self):
        base = "AAATTTCCCGGG"
        mutation = "D:0:A- I:9:-TT I:12:-TTT"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'AATTTCCCTTGGGTTT')
        
    def test_16(self):
        base = "AAATTTCCCGGG"
        mutation = "D:0:A- I:9:-TT D:11:G-"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'AATTTCCCTTGG')
        
    def test_17(self):
        base = "AAATTTCCCGGG"
        mutation = "D:2:ATTT- I:9:-TT"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'AACCCTTGGG')
        
    def test_18(self):
        base = "AAATTTCCCGGG"
        mutation = "D:2:ATTTC- D:10:GG-"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'AACCG')
        
    # Test when front mutant is same base as next base
    def test_19(self):
        base = "AAATTTCCCGGG"
        mutation = "I:0:-A"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'AAAATTTCCCGGG')
    def test_20(self):
        base = "AAATTTCCCGGG"
        mutation = "I:0:-AAA"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'AAAAAATTTCCCGGG')
        
    def test_21(self):
        base = "AAATTTCCCGGG"
        mutation = "I:0:-AAA I:9:-TT"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'AAAAAATTTCCCTTGGG')
        
    def test_22(self):
        base = "AAATTTCCCGGG"
        mutation = "I:0:-AAA I:9:-TT I:12:-TTT"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'AAAAAATTTCCCTTGGGTTT')
        
    def test_23(self):
        base = "AAATTTCCCGGG"
        mutation = "I:0:-AAA I:9:-TT D:11:G-"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'AAAAAATTTCCCTTGG')
        
    def test_24(self):
        base = "AAATTTCCCGGG"
        mutation = "S:3:TTT-AAA I:9:-TT"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'AAAAAACCCTTGGG')
        
    def test_25(self):
        base = "AAATTTCCCGGG"
        mutation = "S:3:TTTC-AAA D:10:GG-"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'AAAAAACCG')
        
    # Test when mutant not in front
    def test_26(self):
        base = "AAATTTCCCGGG"
        mutation = "I:3:-T"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'AAATTTTCCCGGG')
        
    def test_27(self):
        base = "AAATTTCCCGGG"
        mutation = "I:6:-TTCCC"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'AAATTTTTCCCCCCGGG')
        
    def test_28(self):
        base = "AAATTTCCCGGG"
        mutation = "I:6:-A"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'AAATTTACCCGGG')
        
    def test_29(self):
        base = "AAATTTCCCGGG"
        mutation = "I:6:-AA I:9:-AAA"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'AAATTTAACCCAAAGGG')
        
    def test_30(self):
        base = "AAATTTCCCGGG"
        mutation = "S:3:T-A"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'AAAATTCCCGGG')
        
    def test_31(self):
        base = "AAATTTCCCGGG"
        mutation = "S:3:TTT-AAA"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'AAAAAACCCGGG')
    def test_32(self):
        base = "AAATTTCCCGGG"
        mutation = "S:4:T-A"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'AAATATCCCGGG')
        
    def test_33(self):
        base = "AAATTTCCCGGG"
        mutation = "S:6:CCC-AAA"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'AAATTTAAAGGG')
        
    def test_34(self):
        base = "AAATTTCCCGGG"
        mutation = "D:3:T-"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'AAATTCCCGGG')
        
    def test_35(self):
        base = "AAATTTCCCGGG"
        mutation = "D:3:TTT-"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'AAACCCGGG')
        
    def test_36(self):
        base = "AAATTTCCCGGG"
        mutation = "D:5:TC- D:11:G-"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'AAATTCCGG')
        
    def test_37(self):
        base = "AAATTTCCCGGG"
        mutation = "D:3:TTTC- D:11:G-"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'AAACCGG')
        

    # Test tracking
    def test_38(self):
        base = "AAATTTCCCGGG"
        mutation = "I:0:-AAAT D:5:TCC- I:11:-T"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'AAATAAATTCGGTG')
        
    def test_39(self):
        base = "AAATTTCCCGGG"
        mutation = "I:0:-AAAT D:5:T- S:8:C-T D:10:GG-"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'AAATAAATTCCTG')
        
    # Test when mutant behind
    def test_40(self):
        base = "AAATTTCCCGGG"
        mutation = "I:12:-T"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'AAATTTCCCGGGT')
        
    def test_41(self):
        base = "AAATTTCCCGGG"
        mutation = "I:12:-GG"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'AAATTTCCCGGGGG')
        
    def test_42(self):
        base = "AAATTTCCCGGG"
        mutation = "S:11:G-T"
        mutated = mutateSequenceByHash(base, mutation)
        self.assertEqual(mutated, 'AAATTTCCCGGT')
        
    def test_43(self):
        base = "AAATTTCCCGGG"
        mutation = "D:10:GG-"
        mutated = mutateSequenceByHash(base,mutation)
        self.assertEqual(mutated, 'AAATTTCCCG')
        
class MutationExtraction(unittest.TestCase):
    # Test when mutant in front
    def test_1(self):
        base = "AAATTTCCCGGG"
        compare = "TAAATTTCCCGGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'I:0:-T')
        
    def test_2(self):
        base = "AAATTTCCCGGG"
        compare = "TAATTTCCCGGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'S:0:A-T')
        
    def test_3(self):
        base = "AAATTTCCCGGG"
        compare = "AATTTCCCGGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'D:0:A-')
    
    # Test expectation tracking when mutant in front
    # Front mutant is Ins
    def test_4(self):
        base = "AAATTTCCCGGG"
        compare = "TAAATTTCCCTTGGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'I:0:-T I:9:-TT')
        
    def test_5(self):
        base = "AAATTTCCCGGG"
        compare = "TAAATTTCCCTTGGGTTT"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'I:0:-T I:9:-TT I:12:-TTT')
        
    def test_6(self):
        base = "AAATTTCCCGGG"
        compare = "TAAATTTCCCTTGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'I:0:-T I:9:-TT D:11:G-')
        
    def test_7(self):
        base = "AAATTTCCCGGG"
        compare = "TAAACCCTTGGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'I:0:-T D:3:TTT- I:9:-TT')
        
    def test_8(self):
        base = "AAATTTCCCGGG"
        compare = "TAAACCG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'I:0:-T D:3:TTTC- D:10:GG-')
        
    # Test expectation tracking when mutant in front
    # Front mutant is Sub
    def test_9(self):
        base = "AAATTTCCCGGG"
        compare = "TAATTTCCCTTGGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'S:0:A-T I:9:-TT')
        
    def test_10(self):
        base = "AAATTTCCCGGG"
        compare = "TAATTTCCCTTGGGTTT"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'S:0:A-T I:9:-TT I:12:-TTT')
        
    def test_11(self):
        base = "AAATTTCCCGGG"
        compare = "TAATTTCCCTTGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'S:0:A-T I:9:-TT D:11:G-')
        
    def test_12(self):
        base = "AAATTTCCCGGG"
        compare = "TAACCCTTGGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'I:0:-T D:2:ATTT- I:9:-TT')
        
    def test_13(self):
        base = "AAATTTCCCGGG"
        compare = "TAACCG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'I:0:-T D:2:ATTTC- D:10:GG-')

    # Front mutant is Del
    def test_14(self):
        base = "AAATTTCCCGGG"
        compare = "AATTTCCCTTGGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'D:0:A- I:9:-TT')
        
    def test_15(self):
        base = "AAATTTCCCGGG"
        compare = "AATTTCCCTTGGGTTT"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'D:0:A- I:9:-TT I:12:-TTT')
        
    def test_16(self):
        base = "AAATTTCCCGGG"
        compare = "AATTTCCCTTGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'D:0:A- I:9:-TT D:11:G-')
        
    def test_17(self):
        base = "AAATTTCCCGGG"
        compare = "AACCCTTGGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'D:2:ATTT- I:9:-TT')
        
    def test_18(self):
        base = "AAATTTCCCGGG"
        compare = "AACCG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'D:2:ATTTC- D:10:GG-')
        
    # Test when front mutant is same base as next base
    def test_19(self):
        base = "AAATTTCCCGGG"
        compare = "AAAATTTCCCGGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'I:0:-A')
    def test_20(self):
        base = "AAATTTCCCGGG"
        compare = "AAAAAATTTCCCGGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'I:0:-AAA')
        
    def test_21(self):
        base = "AAATTTCCCGGG"
        compare = "AAAAAATTTCCCTTGGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'I:0:-AAA I:9:-TT')
        
    def test_22(self):
        base = "AAATTTCCCGGG"
        compare = "AAAAAATTTCCCTTGGGTTT"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'I:0:-AAA I:9:-TT I:12:-TTT')
        
    def test_23(self):
        base = "AAATTTCCCGGG"
        compare = "AAAAAATTTCCCTTGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'I:0:-AAA I:9:-TT D:11:G-')
        
    def test_24(self):
        base = "AAATTTCCCGGG"
        compare = "AAAAAACCCTTGGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'S:3:TTT-AAA I:9:-TT')
        
    def test_25(self):
        base = "AAATTTCCCGGG"
        compare = "AAAAAACCG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'S:3:TTTC-AAA D:10:GG-')
        
    # Test when mutant not in front
    def test_26(self):
        base = "AAATTTCCCGGG"
        compare = "AAATTTTCCCGGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'I:3:-T')
        
    def test_27(self):
        base = "AAATTTCCCGGG"
        compare = "AAATTTTTCCCCCCGGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'I:6:-TTCCC')
        
    def test_28(self):
        base = "AAATTTCCCGGG"
        compare = "AAATTTACCCGGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'I:6:-A')
        
    def test_29(self):
        base = "AAATTTCCCGGG"
        compare = "AAATTTAACCCAAAGGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'I:6:-AA I:9:-AAA')
        
    def test_30(self):
        base = "AAATTTCCCGGG"
        compare = "AAAATTCCCGGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'S:3:T-A')
        
    def test_31(self):
        base = "AAATTTCCCGGG"
        compare = "AAAAAACCCGGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'S:3:TTT-AAA')
    def test_32(self):
        base = "AAATTTCCCGGG"
        compare = "AAATATCCCGGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'S:4:T-A')
        
    def test_33(self):
        base = "AAATTTCCCGGG"
        compare = "AAATTTAAAGGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'S:6:CCC-AAA')
        
    def test_34(self):
        base = "AAATTTCCCGGG"
        compare = "AAATTCCCGGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'D:3:T-')
        
    def test_35(self):
        base = "AAATTTCCCGGG"
        compare = "AAACCCGGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'D:3:TTT-')
        
    def test_36(self):
        base = "AAATTTCCCGGG"
        compare = "AAATTCCGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'D:5:TC- D:11:G-')
        
    def test_37(self):
        base = "AAATTTCCCGGG"
        compare = "AAACCGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'D:3:TTTC- D:11:G-')
        

    # Test tracking
    def test_38(self):
        base = "AAATTTCCCGGG"
        compare = "AAATAAATTCGGTG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'I:0:-AAAT D:5:TCC- I:11:-T')
        
    def test_39(self):
        base = "AAATTTCCCGGG"
        compare = "AAATAAATTCCTG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'I:0:-AAAT D:5:T- S:8:C-T D:10:GG-')    
        # A better answer
        # self.assertEqual(mutationHash, 'I:4:-AAA S:11:C-T D:12:GG-')    
    
    # Test when mutant behind
    def test_40(self):
        base = "AAATTTCCCGGG"
        compare = "AAATTTCCCGGGT"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'I:12:-T')
        
    def test_41(self):
        base = "AAATTTCCCGGG"
        compare = "AAATTTCCCGGGGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'I:12:-GG')
        
    def test_42(self):
        base = "AAATTTCCCGGG"
        compare = "AAATTTCCCGGT"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'S:11:G-T')
        
    def test_43(self):
        base = "AAATTTCCCGGG"
        compare = "AAATTTCCCG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'D:10:GG-')
        
        
def main():
    unittest.main()

if __name__ == '__main__':
    main()
