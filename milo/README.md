# Requirements
**tqdm**
`pip instal tqdm`
TQDM gives us the progress bar.

**fastcomp**
`pip install fastcomp`
Fast distance check up to distance 2

# How to run

1. Ensure you have the files in the correct folders.
  - `./data/1-raw` should contain matching FASTQ files. 
  - `./references` should contain the Manifest, which are known sequences.

1. Check the config.json for all the configuration options.
  
1. Run `python3 main_ampliconid.py` to convert FASTQ to J3X files
  - Reads fastq files from `1-raw`, outputs j3x to `2-paired`.
  - Unpairable reads are output to `2-paired/discarded`

1. Run `python3 main_mutationid.py` to convert the J3X files to J4X files
  - Reads j3x from `2-paired`, outputs j4x to `3-mutations`.

1. Run `python3 main_mutationstats.py` to gather the stats of the J4X files and create a file to be run by annovar
  - Reads j4x from `3-mutations` and j3x.stats files from `2-paired`, outputs to `4-mutationstats`.

1. Run `python3 main_annovar.py` to run against annovar to check for existing mutations

# File format documentation

[Same as the algo documentation](https://docs.google.com/document/d/1_uWV8ExxDhpnAHwQIGdE2CcQR7scXawhVlBf9aegF8Q/edit?usp=sharing)