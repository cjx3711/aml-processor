# Requirements
**tqdm**
`pip instal tqdm`
TQDM gives us the progress bar.

**fastcomp**
`pip install fastcomp`
Fast distance check up to distance 2

# How to run

1. Ensure you have the files in the correct folders.
  - `./data/Raw` should contain matching FASTQ files. 
  - `./references` should contain the matching 
  
1. Run `python3 main-ampliconid.py` to convert FASTQ to J3X files

1. Run `python3 main-mutationid.py` to convert the J3X files to J4X files

1. Run `python3 main-mutationstats.py` to gather the stats of the J4X files and create a file to be run by annovar

# Check against annovar

1. Install annovar into a directory called `annovar` next to `milo`
1. Run annovar, which will output a file `myanno.hg19_multianno.csv`

**Running from the `annovar` dir**
```
perl table_annovar.pl ../milo/4-mutationstats/annovarStats.csv humandb/ -buildver hg19 -out myanno -remove -protocol refGene,avsnp147,clinvar_20170130,cosmic70 -operation gx,f,f,f -nastring . -csvout -polish -xref example/gene_xref.txt
```

**Running from the `milo` dir**
Next run the script to combine both files
```
python main_annovar.py data/4-mutationstats/annovarStats.csv ../annovar/myanno.hg19_multianno.csv
```

# File format documentation

[Same as the algo documentation](https://docs.google.com/document/d/1_uWV8ExxDhpnAHwQIGdE2CcQR7scXawhVlBf9aegF8Q/edit?usp=sharing)