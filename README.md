# MILo: Mutant Indel Locator

[![Build Status](https://travis-ci.org/j3x1/aml-processor.svg?branch=master)](https://travis-ci.org/j3x1/aml-processor)

MILo is a FastQ Aligner and Variant Caller algorithm for the AML Panel.

**Stack**
- Python 3.6

[Algo documentation](https://docs.google.com/document/d/1_uWV8ExxDhpnAHwQIGdE2CcQR7scXawhVlBf9aegF8Q/edit?usp=sharing).

## Jargon
```
Unpaired left read:
ATCGATCGTTATCG

Unpaired right read:
CGATCGATAACGAT
       ↓
RC'ed right read:
ATCGTTATCGATCG

(Proc) Align pairs:
ATCGATCGTTATCG
    ATCGTTATCGATCG
       ↓
Paired read:
ATCGATCGTTATCGATCG

Signed paired reads:
    Negative paired read:
    - CGATCGATAACGATCGAT
              ↓
    Positive paired read:
    + ATCGATCGTTATCGATCG

Unsigned paired reads:
    (Proc) _Tile paired reads_:
    ATCGATCGTTATCGATCG
              ATCGATCGTTATCGATCG
                        ATCGATCGTTATCGATCG
                    ↓
    _Tiled reads_:
    ATCGATCGTTATCGATCGTTATCGATCGTTATCGATCG
```
1. Pair reads, 3 possibilities:
   1) Aligns.
      a) Subs/1bp indels
      b) > 1bp insertions
      c) > 1bp deletions
   2) Alignment failure due to overlap exceeding line (large insertion): mark header with ¦, insert … between reads.
   3) Alignment failure due to low quality substrings: mark header with ‽, insert space between reads.

   WIP Algo:
   Reverse complement right read.
   Check distance between left and right reads. If distance ≤ 5, then a.
   Else, get high quality substrings (break up if > 10bp) from left read and insert into list. Search for the first reasonably long substring from both ends of the list find them in the right read. If can't find, try again with the next reasonably long substring.

2. Align Reads:
   Checks the overlap between left and right reads for differences, and picks the nucleotide with a higher quality score. Combines both reads into 1. Delete probe ends if necessary.
      Header left with the coordinates, and ¦ or ‽.
      Reformats quality scores as follows:
      100% > confidence ≥ 99.9%:       .
      99.9% > confidence ≥ 99%:        ,
      99% > confidence ≥ 90%:          ?
      0% > confidence ≥ 90%:           ×

   Match to amplicon, reverse complement if direction is -

   Match to combined sequenced region
