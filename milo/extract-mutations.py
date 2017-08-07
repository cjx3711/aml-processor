from j4xUtils import *

base = "AAATTTCCCGGG" 
seq1 = "AAATTTCCUGGG" # sub U
seq2 = "AAATTTCCGGGG" # sub G
seq3 = "AAAAATTTCCCGGGG" # ins AA and G
seq4 = "AAATTTTTTCCCGGG" # ins TTT
seq5 = "AATTCCCGGG" # del AT
seq6 = "AAAUUTTTTACCCGH"
seq7 = "TTTGGG" # del AAA del CCC
seq8 = "AAATTTCCCGGC"
seq9 = "AAATTTCCATATGGG" 


mutationID(base, seq1)
mutationID(base, seq2)
mutationID(base, seq3)
mutationID(base, seq4)
mutationID(base, seq5)
mutationID(base, seq6)
mutationID(base, seq7)
mutationID(base, seq8)
mutationID(base, seq9)