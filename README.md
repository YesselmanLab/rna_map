


# RNA MAP

[![Docker Linux Build](https://github.com/YesselmanLab/rna_map/actions/workflows/docker_linux_build.yml/badge.svg)](https://github.com/YesselmanLab/rna_map/actions/workflows/docker_linux_build.yml)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![linting: flake8](https://img.shields.io/badge/linting-flake8-greenyellow)](https://github.com/PyCQA/flake8)
[![PYPI package](https://badge.fury.io/py/rna-map.png)](http://badge.fury.io/py/rna-map)

A open-source tool for rapid analysis of RNA mutational profiling (MaP) experiments.
This tool was inspired by the DREEM algorithm developed by the Rouskin Lab (https://www.rouskinlab.com/). Please cite this work (https://doi.org/10.1093/nar/gkac435).

The MaP analysis web tool provides a simple platform for analyzing DMS-reactivity of an RNA. The user input is a raw sequencing file (.fastq) generated from a DMS-MaPseq experiment, and a sequence of the RNA of interest (.fasta). The DREEM algorithm performs sequence alignment using bowtie-2 and outputs the mismatch rate per nucleotide.

## Software requirements

- python 3.8 or greater
- bowtie2 (2.2.9) - https://github.com/BenLangmead/bowtie2/releases/download/v2.2.9/
- trim_galore (0.6.6) - https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz
- fastqc (0.11.9) - https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
- cutadapt (1.18) - https://github.com/marcelm/cutadapt/archive/refs/tags/v1.18.zip
- conda (optional) - https://docs.anaconda.com/anaconda/install/
- docker (optional) - https://docs.docker.com/get-docker/

If you are trying the software for the first time highly recommended to use the docker image.

## How to install
Highly recommended to use conda to manage your python environment.

```
conda create -n rna-map python=3.8
pip install rna-map
```

### with docker 
```shell
# on linux and intel mac
git clone https://github.com/YesselmanLab/rna_map
cd rna_map
pip install .

docker build -t rna-map -f docker/Dockerfile .

# on mac with apple silicon / or other arm64 platforms
docker build -t rna-map --platform linux/amd64 -f docker/Dockerfile .

```

## How to use 

### basic usage

After installed there will be a new command line tool called `rna-map` available.

```shell
# run in single end mode
rna-map -fa <fasta file> -fq1 <fastq file>

# run in paired end mode
rna-map -fa <fasta file> -fq1 <fastq file> -fq2 <fastq file>

# supply a csv with dot bracket structures. These will apppear in the 
# results in plots and pickle file 
rna-map -fa <fasta file> -fq1 <fastq file> --dot-bracket <csv file>
```

### running with docker 
`--docker` flag will run the docker image. if you have run docker build first
```shell
# run in single end mode
# note this will only work if you built the image with docker command above
rna-map -fa <fasta file> -fq1 <fastq file> --docker

# TODO check is necessary? 
# run on apple silicon on / or other arm64 platforms 
rna-map -fa <fasta file> -fq1 <fastq file> --docker --docker-platform linux/amd64
```

### working with large sets of RNAs



### see a full list of arguments below

```
rna-map --help
Usage: rna-map [OPTIONS]

  rapid analysis of RNA mutational profiling (MaP) experiments.

Main arguments:
  These are the main arguments for the command line interface
  -fa, --fasta PATH              The fasta file containing the reference
                                 sequences  [required]
  -fq1, --fastq1 PATH            The fastq file containing the single end reads
                                 or the first pair of paired end reads
                                 [required]
  -fq2, --fastq2 TEXT            The fastq file containing the second pair of
                                 paired end reads
  --dot-bracket TEXT             The directory containing the input files
  -pf, --param-file TEXT         A yml formatted file to specify parameters, see
                                 rna_map/resources/default.yml for an example
  -pp, --param-preset TEXT       run a set of parameters for specific uses like
                                 'barcoded-libraries'

Mapping options:
  These are the options for pre processing of fastq files and alignment to
  reference sequences
  --skip-fastqc                  do not run fastqc for quality control of
                                 sequence data
  --skip-trim-galore             do not run trim galore for quality control of
                                 sequence data
  --tg-q-cutoff INTEGER          the quality cutoff for trim galore
  --bt2-alignment-args TEXT      the arguments to pass to bowtie2 for alignment
                                 seperated by commas
  --save-unaligned               the path to save unaligned reads to

Bit vector options:
  These are the options for the bit vector step
  --skip-bit-vector              do not run the bit vector step
  --summary-output-only          do not generate bit vector files or plots
                                 recommended when there are thousands of
                                 reference sequences
  --plot-sequence                plot sequence and structure is supplied under
                                 the population average plots
  --map-score-cutoff INTEGER     reject any bit vector where the mapping score
                                 for bowtie2 alignment is less than this value
  --qscore-cutoff INTEGER        quality score of read nucleotide, sets to
                                 ambigious if under this val
  --mutation-count-cutoff INTEGER
                                 maximum number of mutations allowed in a bit
                                 vector will be discarded if higher
  --percent-length-cutoff FLOAT  minium percent of the length of the reference
                                 sequence allowed in a bit vector will be
                                 discarded if lower
  --min-mut-distance INTEGER     minimum distance between mutations in a bit
                                 vector will be discarded if lower

Docker options:
  These are the options for running the command line interface in a docker
  container
  --docker                       Run the program in a docker container
  --docker-image TEXT            The docker image to use
  --docker-platform TEXT         The platform to use for the docker image

Misc options:
  These are the options for the misc stage
  --overwrite                    overwrite the output directory if it exists
  --restore-org-behavior         restore the original behavior of the rna_map
  --stricter-bv-constraints      use stricter constraints for bit vector
                                 generation, use at your own risk!
  --debug                        enable debug mode

Other options:
  --help                         Show this message and exit.

```

### running paired end reads

```shell
 rna-map -fa test/resources/case_1/test.fasta -fq1 test/resources/case_unit/test_mate1.fastq -fq2 test/resources/case_unit/test_mate2.fastq 
```

```shell
rna_map.CLI - INFO -
88888888ba   888b      88         db             88b           d88         db         88888888ba
88      "8b  8888b     88        d88b            888b         d888        d88b        88      "8b
88      ,8P  88 `8b    88       d8'`8b           88`8b       d8'88       d8'`8b       88      ,8P
88aaaaaa8P'  88  `8b   88      d8'  `8b          88 `8b     d8' 88      d8'  `8b      88aaaaaa8P'
88""""88'    88   `8b  88     d8YaaaaY8b         88  `8b   d8'  88     d8YaaaaY8b     88""""""'
88    `8b    88    `8b 88    d8""""""""8b        88   `8b d8'   88    d8""""""""8b    88
88     `8b   88     `8888   d8'        `8b       88    `888'    88   d8'        `8b   88
88      `8b  88      `888  d8'          `8b      88     `8'     88  d8'          `8b  88

rna_map.CLI - INFO - ran at commandline as:
rna_map.CLI - INFO - /Users/jyesselm/miniconda3/envs/py3/bin/rna-map -fa test/resources/case_1/test.fasta -fq1 test/resources/case_unit/test_mate1.fastq -fq2 test/resources/case_unit/test_mate2.fastq
rna_map.RUN - INFO - fasta file: test/resources/case_1/test.fasta exists
rna_map.RUN - INFO - found 1 valid reference sequences in test/resources/case_1/test.fasta
rna_map.RUN - INFO - fastq1 file: test/resources/case_unit/test_mate1.fastq exists
rna_map.RUN - INFO - fastq2 file: test/resources/case_unit/test_mate2.fastq exists
rna_map.RUN - INFO - two fastq files supplied, thus assuming paired reads
rna_map.MAPPING - INFO - bowtie2 2.4.5 detected!
rna_map.MAPPING - INFO - fastqc v0.11.9 detected!
rna_map.MAPPING - INFO - trim_galore 0.6.6 detected!
rna_map.MAPPING - INFO - cutapt 1.18 detected!
rna_map.MAPPING - INFO - building directory structure
rna_map.MAPPING - INFO - bowtie2 2.4.5 detected!
rna_map.MAPPING - INFO - fastqc v0.11.9 detected!
rna_map.MAPPING - INFO - trim_galore 0.6.6 detected!
rna_map.MAPPING - INFO - cutapt 1.18 detected!
rna_map.EXTERNAL_CMD - INFO - running fastqc
rna_map.EXTERNAL_CMD - INFO - fastqc ran without errors
rna_map.EXTERNAL_CMD - INFO - running trim_galore
rna_map.EXTERNAL_CMD - INFO - trim_galore ran without errors
rna_map.EXTERNAL_CMD - INFO - running bowtie2-build
rna_map.EXTERNAL_CMD - INFO - bowtie2-build ran without errors
rna_map.EXTERNAL_CMD - INFO - running bowtie2 alignment
rna_map.EXTERNAL_CMD - INFO - bowtie2 alignment ran without errors
rna_map.EXTERNAL_CMD - INFO - results for bowtie alignment:
25 reads; of these:
  25 (100.00%) were paired; of these:
    1 (4.00%) aligned concordantly 0 times
    24 (96.00%) aligned concordantly exactly 1 time
    0 (0.00%) aligned concordantly >1 times
96.00% overall alignment rate
rna_map.MAPPING - INFO - finished mapping!
rna_map.BIT_VECTOR - INFO - starting bitvector generation
rna_map.BIT_VECTOR - INFO - REMOVED READS:
| name          |   low_mapq |
|---------------|------------|
| mttr-6-alt-h3 |          0 |

rna_map.BIT_VECTOR - INFO - MUTATION SUMMARY:
| name          |   reads |   aligned |   no_mut |   1_mut |   2_mut |   3_mut |   3plus_mut |   sn |
|---------------|---------|-----------|----------|---------|---------|---------|-------------|------|
| mttr-6-alt-h3 |      24 |       100 |       50 |   33.33 |    12.5 |    4.17 |           0 | 4.91 |
```



## TODO
- [ ] 
- [ ] Add mac build to github actions
