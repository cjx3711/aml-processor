import unittest

from genomicsUtils import *
from pairedFastqProc import *



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


class PairedFASTQAligher(unittest.TestCase):


    def test_no_overlap(self):
        left = "A" * 151
        right = "C" * 151
        leftQuality = "K" * len(left)
        rightQuality = "K" * len(right)

        mergedSequence, qualityScores, noOfCol = mergeUnpaired(left, right, leftQuality, rightQuality)
        expected = left + " " + right

        # No overlap because too small
        self.assertEqual(mergedSequence, expected)
        self.assertEqual(noOfCol, "N/A")

    def test_small_overlap(self):
        overlap = 5
        sides = 146
        left = "A" * sides + "C" * overlap
        right = "C" * overlap + "G" * sides
        leftQuality = "K" * len(left)
        rightQuality = "K" * len(right)

        mergedSequence, qualityScores, noOfCol = mergeUnpaired(left, right, leftQuality, rightQuality)

        expected = left + " " + right

        # No overlap because too small
        self.assertEqual(mergedSequence, expected)
        self.assertEqual(noOfCol, "N/A")

    def test_normal_overlap(self):
        overlap = 11
        sides = 140
        left = "A" * sides + "C" * overlap
        right = "C" * overlap + "G" * sides
        leftQuality = "K" * len(left)
        rightQuality = "K" * len(right)

        mergedSequence, qualityScores, noOfCol = mergeUnpaired(left, right, leftQuality, rightQuality)

        expected = "A" * sides + "C" * overlap + "G" * sides

        self.assertEqual(mergedSequence, expected)

    def test_large_overlap(self):
        overlap = 51
        sides = 100
        left = "A" * sides + "C" * overlap
        right = "C" * overlap + "G" * sides
        leftQuality = "K" * len(left)
        rightQuality = "K" * len(right)

        mergedSequence, qualityScores, noOfCol = mergeUnpaired(left, right, leftQuality, rightQuality)

        expected = "A" * sides + "C" * overlap + "G" * sides

        # vulnerable to homopolymer sequences. It should give the wrong value if the overlap is too big
        self.assertNotEqual(mergedSequence, expected)

    def test_actual_sequence(self):
        left = "CAAGTTGATGGGGGCAAAACCAAAATAATTTTCATTTTTAATATACCACACAACACATCTATCTACAAATGCTTTACATTAAATTTACCTGCCAACTGTTTAGCCTGGCTTGCGTTTTCAGTTTGTTTTGTACGTGATGGGGCTGACTTTT"
        right = reverseComplement("AGTCAGGATGTTAGCAGAGCCAGTCAAGACTTGCCGACAAAGGAAACTAGAAGCCAAGAAAGCTGCAGCTGAAAAGCTTTCCTCCCTGGAGAACAGCTCAAATAAAAATGAAAAGGAAAAGTCAGCCCCATCACGTACAAAACAAACTGAA")
        lQuality = "CCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
        rQuality = "CCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGDGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGDFG"[::-1]
        mergedSequence, qualityScores, noOfCol = mergeUnpaired(left, right, lQuality, rQuality)

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

        mergedSequence, qualityScores, noOfCol = mergeUnpaired(left, right, lQuality, rQuality)

        expectedMerged = "TCTTTCTGCCTCATGCTCTCTCCAACAGGCTTGCAGCCAATTTACTGGAGCAGGGATGACGTAGCCCAGTGGCTCATGTGGGCTGAAAATGAGTTTTCTTTAAGGCCAATTGACAGCAACACATTTGAAATGAATGGCAAAGCTCTCCTGCTGCTGACCAAAGAGGACTTTCGCTATCGATCTCCTCATTCAGGTGAGAGTCTGGACTCTTGGCATATGCCCAACTTGGAAAATCTCTTAGTTAGTGGTTGGTCTTTAACACCCC"
        expecedQuality = ".............,.,...................,..,.........,.....,,.,..,..,,.....,.....?..........,.,.,.,,,......................,,.?..............,.......,.,,...,?.........,?,?.,,...,.......,...?,.?,?,,??.?,....,,.?..?,....,...?......,???....?.......,.........,??.,,....,..,?"
        expectedCol = '0'
        self.assertEqual(mergedSequence, expectedMerged)
        self.assertEqual(qualityScores, expecedQuality)
        self.assertEqual(noOfCol, expectedCol)

    def test_actual_sequence_delete(self):
        left = "TGCAGTCCCAGCCCACAGCCCCCCTCCTCCCTCAGACTCAGGAGTCCATAGCGAATTTCGACGATCGTTGCATTAACTCGCGAACAAACGGATCTCGTATGCCGTCTTCTGCTTGAAAACAAAAACAATTTCAATACGTCTGAAAATGTTA"
        right = reverseComplement("TATGGACTCCTGAGTCTGAGGGAGGAGGGGGGCTGTGGGCTGGGACTGCAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAAGGTTCAGTGTAGATCTCGGTGGGCGCCGGATCATTAATAAAAAAAAAACAGAAATTCACGGTGAGC")
        print(right)
        lQuality = "CCCCCGEEFGGGGGGGGEGCCEFGCEGGFCFGGGGDG@EFGFGGGAFGGGEGGGGCDEEGFGGEEEFEG@CFEDEEGCFGCEFE@FEGFGG=4+,EFF7D=FGG=CEFB??,CFGFGG<,?,,4B+:,,,8<5,,C,B:,,+5,:,:,@B,"
        rQuality = "CCCCCFGGGGGGFGGFGFFGGGGGG;CEDEGGGGGGFGGGGGDG<AEFG?EBFFGGG:@CEGCF>FECFFFGGE<BFF,BCE,EABFB:>>FDFF9FF,=3@F<>:C*3<*5<****8,<;<,@,8;:=E*1*1*++2,27,,,***4*2+"[::-1]

        mergedSequence, qualityScores, noOfCol = mergeUnpaired(left, right, lQuality, rQuality)

        expectedMerged = "TGCAGTCCCAGCCCACAGCCCCCCTCCTCCCTCAGACTCAGGAGTCCATA"
        expecedQuality = ".................................................."
        expectedCol = '0'
        self.assertEqual(mergedSequence, expectedMerged)
        self.assertEqual(qualityScores, expecedQuality)
        self.assertEqual(noOfCol, expectedCol)


def main():
    unittest.main()

if __name__ == '__main__':
    main()
