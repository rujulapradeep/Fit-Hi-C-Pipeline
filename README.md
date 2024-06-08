# Fit-Hi-C Pipeline
 A python pipeline designed to generate a comprehensive list of intra and inter significant chromosomal interactions for Hi-C data.

 # Pre-requirements
 - Python
 - JuicerTools: https://github.com/aidenlab/JuicerTools
 - fithic: https://github.com/ay-lab/fithic
 - gzip
 - shutil
 - glob
 - argparse
   
# Instructions
In order to execute the pipeline, you will need to first install the requirements above. Additionally, you will need a list of chromosome sizes.

## Sample Data Files
### Example Hi-C Sample File
```
[path_to_hic_sample] [sample_name]
path_to_hic_sample1 sample1
path_to_hic_sample2 sample2
path_to_hic_sample3 sample3
path_to_hic_sample4 sample4
path_to_hic_sample5 sample5
```
#### *Important*
- two columns
- separated by space
- no header

### Chromosome Sizes List

```
chr1 ####
chr2 ####
chr3 ####
chr4 ####
chr5 ####
```
#### *Important*
- two columns
- separated by space
- no header

# fithic.py
```
usage: fithic.py [-h] <hic_samples_fil> <output_dir> <chromsizes_file> <juicer_jar_path> <fithic_path> <resolution>

Run the Fit-Hi-C pipeline for Hi-C data.

positional arguments:
  hic_samples_file  Path to the text file containing Hi-C sample paths and names
  output_dir        Directory to save output files
  chromsizes_file   Path to the chromosome sizes file
  juicer_jar_path   Path to Juicer Tools .jar file
  fithic_path       Path to the Fit-Hi-C directory
  resolution        Resolution for the Fit-Hi-C analysis

optional arguments:
  -h, --help        show this help message and exit

```
## Output
You will find a comprehensive list of both intra- and inter-chromosomal significant interactions within the output directory, specifically under the Significant_Interactions directory, organized by each Hi-C sample.
