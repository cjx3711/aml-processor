import unittest

from genomicsUtils import *
from j4xUtils import *
from AmpliconMatcherProbabilistic import *
from AmpliconMatcherHashSweep import *
from ReadPairer import *

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

class AmpliconMatcherProbabilisticTests(unittest.TestCase):
    def test_file_read(self):
        # Count the amplicons that were read
        ampMat = AmpliconMatcherProbabilistic("references/Manifest.csv")
        self.assertEqual(len(ampMat.refSeqs), 571)
        self.assertEqual(len(ampMat.ampDict1), 571)
        self.assertEqual(len(ampMat.ampDict2), 1713)
        self.assertEqual(len(ampMat.refSeqsTrunc), 571)

class AmpliconMatcherHashSweepTests(unittest.TestCase):
    def test_file_read(self):
        # Count the amplicons that were read
        ampMat = AmpliconMatcherHashSweep("references/Manifest.csv")
        self.assertEqual(ampMat.referenceCount, 571)

    def test_simple_data(self):
        sample_data = [
            'AAAATTTT',
            'CCCCGGGG'
        ]
        ampMat = AmpliconMatcherHashSweep(sample_data, 4, 4)
        self.assertEqual(ampMat.referenceCount, 2)

        # Check if the hash structure is correct
        self.assertTrue('AAAA' in ampMat.ampliconRefs)
        self.assertTrue('TTTT' in ampMat.ampliconRefs)
        self.assertTrue('CCCC' in ampMat.ampliconRefs)
        self.assertTrue('GGGG' in ampMat.ampliconRefs)
        self.assertTrue((0, 0) in ampMat.ampliconRefs['AAAA'])
        self.assertTrue((0, 4) in ampMat.ampliconRefs['TTTT'])
        self.assertTrue((1, 0) in ampMat.ampliconRefs['CCCC'])
        self.assertTrue((1, 4) in ampMat.ampliconRefs['GGGG'])

    def test_overlapping_data(self):
        sample_data = [
            'AAAATTTT',
            'TTTTGGGG'
        ]
        ampMat = AmpliconMatcherHashSweep(sample_data, 4, 4)
        self.assertEqual(ampMat.referenceCount, 2)

        # Check if the hash structure is correct
        self.assertTrue('AAAA' in ampMat.ampliconRefs)
        self.assertTrue('TTTT' in ampMat.ampliconRefs)
        self.assertFalse('CCCC' in ampMat.ampliconRefs)
        self.assertTrue('GGGG' in ampMat.ampliconRefs)
        self.assertEqual(len(ampMat.ampliconRefs['TTTT']), 2)
        self.assertTrue((0, 0) in ampMat.ampliconRefs['AAAA'])
        self.assertTrue((0, 4) in ampMat.ampliconRefs['TTTT'])
        self.assertTrue((1, 0) in ampMat.ampliconRefs['TTTT'])
        self.assertTrue((1, 4) in ampMat.ampliconRefs['GGGG'])

    def test_skip_data(self):
        sample_data = [
            'AAAATTTTCCCCGGGG'
        ]
        ampMat = AmpliconMatcherHashSweep(sample_data, 4, 5)
        self.assertEqual(ampMat.referenceCount, 1)

        # Check if the hash structure is correct
        self.assertTrue('AAAA' in ampMat.ampliconRefs)
        self.assertFalse('TTTT' in ampMat.ampliconRefs)
        self.assertFalse('CCCC' in ampMat.ampliconRefs)
        self.assertFalse('GGGG' in ampMat.ampliconRefs)
        self.assertTrue('TTTC' in ampMat.ampliconRefs)
        self.assertTrue('CCGG' in ampMat.ampliconRefs)
        self.assertTrue((0, 0) in ampMat.ampliconRefs['AAAA'])
        self.assertTrue((0, 5) in ampMat.ampliconRefs['TTTC'])
        self.assertTrue((0, 10) in ampMat.ampliconRefs['CCGG'])



class PairedFASTQAligner(unittest.TestCase):

    def setUp(self):
        self.readPairer = ReadPairer()

    def test_scoring(self):
        left = "AAAAAAAAAAAAA"
        right = "AAAAAAAAAAAAA"
        overlapPairs = tuple(zip(left, right))
        score = self.readPairer.calcScore(overlapPairs, len(left))
        self.assertEqual(score, 13)

    def test_score_ignore(self):
        left = "AAAAAAAAAAAAA"
        right = "CCCCCCCCCCCCC"
        overlapPairs = tuple(zip(left, right))
        score = self.readPairer.calcScore(overlapPairs, len(left))
        self.assertEqual(score, -5)

    def test_score_mismatch(self):
        left = "AAAAAAAAAAAAA"
        right = "AAAAAAAAAAAAC"
        overlapPairs = tuple(zip(left, right))
        score = self.readPairer.calcScore(overlapPairs, len(left))
        self.assertEqual(score, 7)

    def test_no_overlap(self):
        left = "A" * 151
        right = "C" * 151
        leftQuality = "K" * len(left)
        rightQuality = "K" * len(right)

        mergedSequence, qualityScores, noOfCol = self.readPairer.mergeUnpaired(left, right, leftQuality, rightQuality)
        expected = left + " " + right

        # No overlap because too small
        self.assertEqual(mergedSequence, expected)
        self.assertEqual(noOfCol, "?")

    def test_small_overlap(self):
        overlap = 5
        sides = 146
        left = "A" * sides + "C" * overlap
        right = "C" * overlap + "G" * sides
        leftQuality = "K" * len(left)
        rightQuality = "K" * len(right)

        mergedSequence, qualityScores, noOfCol = self.readPairer.mergeUnpaired(left, right, leftQuality, rightQuality)

        expected = left + " " + right

        # No overlap because too small
        self.assertEqual(mergedSequence, expected)
        self.assertEqual(noOfCol, "?")

    def test_normal_overlap(self):
        overlap = 11
        sides = 140
        left = "A" * sides + "C" * overlap
        right = "C" * overlap + "G" * sides
        leftQuality = "K" * len(left)
        rightQuality = "K" * len(right)

        mergedSequence, qualityScores, noOfCol = self.readPairer.mergeUnpaired(left, right, leftQuality, rightQuality)

        expected = "A" * sides + "C" * overlap + "G" * sides

        self.assertEqual(mergedSequence, expected)

    def test_large_overlap(self):
        overlap = 51
        sides = 100
        left = "A" * sides + "C" * overlap
        right = "C" * overlap + "G" * sides
        leftQuality = "K" * len(left)
        rightQuality = "K" * len(right)

        mergedSequence, qualityScores, noOfCol = self.readPairer.mergeUnpaired(left, right, leftQuality, rightQuality)

        expected = "A" * sides + "C" * overlap + "G" * sides

        # vulnerable to homopolymer sequences. It should give the wrong value if the overlap is too big
        self.assertNotEqual(mergedSequence, expected)

    def test_actual_sequence(self):
        left = "CAAGTTGATGGGGGCAAAACCAAAATAATTTTCATTTTTAATATACCACACAACACATCTATCTACAAATGCTTTACATTAAATTTACCTGCCAACTGTTTAGCCTGGCTTGCGTTTTCAGTTTGTTTTGTACGTGATGGGGCTGACTTTT"
        right = reverseComplement("AGTCAGGATGTTAGCAGAGCCAGTCAAGACTTGCCGACAAAGGAAACTAGAAGCCAAGAAAGCTGCAGCTGAAAAGCTTTCCTCCCTGGAGAACAGCTCAAATAAAAATGAAAAGGAAAAGTCAGCCCCATCACGTACAAAACAAACTGAA")
        lQuality = "CCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
        rQuality = "CCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGDGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGDFG"[::-1]
        mergedSequence, qualityScores, noOfCol = self.readPairer.mergeUnpaired(left, right, lQuality, rQuality)

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

        mergedSequence, qualityScores, noOfCol = self.readPairer.mergeUnpaired(left, right, lQuality, rQuality)

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

        mergedSequence, qualityScores, noOfCol = self.readPairer.mergeUnpaired(left, right, lQuality, rQuality)

        expectedMerged = "TGCAGTCCCAGCCCACAGCCCCCCTCCTCCCTCAGACTCAGGAGTCCATA"
        expecedQuality = ".................................................."
        expectedCol = '0'
        self.assertEqual(mergedSequence, expectedMerged)
        self.assertEqual(qualityScores, expecedQuality)
        self.assertEqual(noOfCol, expectedCol)


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
        self.assertEqual(mutationHash, 'I:0:-T I:10:-TT')
        
    def test_5(self):
        base = "AAATTTCCCGGG"
        compare = "TAAATTTCCCTTGGGTTT"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'I:0:-T I:10:-TT I:15:-TTT')
        
    def test_6(self):
        base = "AAATTTCCCGGG"
        compare = "TAAATTTCCCTTGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'I:0:-T I:10:-TT D:14:G-')
        
    def test_7(self):
        base = "AAATTTCCCGGG"
        compare = "TAAACCCTTGGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'I:0:-T D:4:TTT- I:7:-TT')
        
    def test_8(self):
        base = "AAATTTCCCGGG"
        compare = "TAAACCG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'I:0:-T D:4:TTTC- D:7:GG-')
        
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
        self.assertEqual(mutationHash, 'S:0:A-T I:9:-TT I:14:-TTT')
        
    def test_11(self):
        base = "AAATTTCCCGGG"
        compare = "TAATTTCCCTTGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'S:0:A-T I:9:-TT D:13:G-')
        
    def test_12(self):
        base = "AAATTTCCCGGG"
        compare = "TAACCCTTGGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'I:0:-T D:3:ATTT- I:6:-TT')
        
    def test_13(self):
        base = "AAATTTCCCGGG"
        compare = "TAACCG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'I:0:-T D:3:ATTTC- D:6:GG-')

    # Front mutant is Del
    def test_14(self):
        base = "AAATTTCCCGGG"
        compare = "AATTTCCCTTGGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'D:0:A- I:8:-TT')
        
    def test_15(self):
        base = "AAATTTCCCGGG"
        compare = "AATTTCCCTTGGGTTT"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'D:0:A- I:8:-TT I:13:-TTT')
        
    def test_16(self):
        base = "AAATTTCCCGGG"
        compare = "AATTTCCCTTGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'D:0:A- I:8:-TT D:12:G-')
        
    def test_17(self):
        base = "AAATTTCCCGGG"
        compare = "AACCCTTGGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'D:2:ATTT- I:5:-TT')
        
    def test_18(self):
        base = "AAATTTCCCGGG"
        compare = "AACCG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'D:2:ATTTC- D:5:GG-')
        
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
        self.assertEqual(mutationHash, 'I:0:-AAA I:12:-TT')
        
    def test_22(self):
        base = "AAATTTCCCGGG"
        compare = "AAAAAATTTCCCTTGGGTTT"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'I:0:-AAA I:12:-TT I:17:-TTT')
        
    def test_23(self):
        base = "AAATTTCCCGGG"
        compare = "AAAAAATTTCCCTTGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'I:0:-AAA I:12:-TT D:16:G-')
        
    def test_24(self):
        base = "AAATTTCCCGGG"
        compare = "AAAAAACCCTTGGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'S:3:TTT-AAA I:9:-TT')
        
    def test_25(self):
        base = "AAATTTCCCGGG"
        compare = "AAAAAACCG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'S:3:TTTC-AAA D:9:GG-')
        
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
        self.assertEqual(mutationHash, 'I:6:-AA I:11:-AAA')
        
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
        self.assertEqual(mutationHash, 'D:5:TC- D:9:G-')
        
    def test_37(self):
        base = "AAATTTCCCGGG"
        compare = "AAACCGG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'D:3:TTTC- D:7:G-')
        

    # Test tracking
    def test_38(self):
        base = "AAATTTCCCGGG"
        compare = "AAATAAATTCGGTG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'I:0:-AAAT D:9:TCC- I:12:-T')
        
    def test_39(self):
        base = "AAATTTCCCGGG"
        compare = "AAATAAATTCCTG"
        mutationHash = mutationArrayToHash(mutationID(base,compare))
        self.assertEqual(mutationHash, 'I:0:-AAAT D:9:T- S:11:C-T D:13:GG-')    
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
