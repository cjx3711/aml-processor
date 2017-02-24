# FastQ Aligner and Variant Caller for AML Panel
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

```
Checks the overlap between left and right reads for differences, and picks the nucleotide with a higher quality score. Combines both reads into 1. Delete probe ends if necessary.
	Header left with the coordinates, and ¦ or ‽.
	Reformats quality scores as follows:
	100% > confidence ≥ 99.9%:       .
	99.9% > confidence ≥ 99%:        ,
	99% > confidence ≥ 90%:          ?
	0% > confidence ≥ 90%:           ×

Match to amplicon, reverse complement if direction is -

Match to combined sequenced region
```