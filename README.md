[![Build Status](https://app.travis-ci.com/jyesselm/dreem.svg?branch=main)](https://app.travis-ci.com/jyesselm/dreem)


# DREEM 

DREEM processes DMS next-generation sequencing data to produce mutational profiles that relate to DMS modification rates written by Silvi Rouskin and the Rouskin lab (https://www.rouskinlab.com/)

## Install 

The easiest way to install is via docker! 

If you are not familiar with docker https://docs.docker.com/get-docker/

```shell
# get the code and install the python module
wget https://github.com/jyesselm/dreem/archive/refs/tags/0.2.0.zip      && \
    unzip 0.2.0.zip                                                     && \ 
    mv dreem-0.2.0 dreem                                                && \
    cd dreem                                                            && \
    sudo python3 setup.py install                                       && \
    pip install .                                                       && \
    pip install -r requirements.txt

# to setup the docker image you can either pull
docker pull jyesselm/dreem:latest
# or build docker image from scratch 
docker build -t dreem -f docker/Dockerfile .

```

If you want to run natively DREEM requires:

```she
bowtie2=v2.2.9 (https://github.com/BenLangmead/bowtie2/releases/download/v2.2.9/)
fastqc (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
cutapt=v1.18 (https://github.com/marcelm/cutadapt/archive/refs/tags/v1.18.zip)
TrimGalore=v0.6.6 (https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz)
```

I also highly recommend using anaconda for installing all required python packages which are all in requirements.txt

```shell
colorlog
click
click-option-group
biopython
plotly
matplotlib
pandas
tabulate
pyyaml
```

## Run 

### running locally 

main excutable is `dreem` after installing locally it should be available in $PATH 

```shell
Usage: run.py [OPTIONS]

  DREEM processes DMS next generation sequencing data to produce mutational
  profiles that relate to DMS modification rates written by Silvi Rouskin and
  the Rouskin lab (https://www.rouskinlab.com/)

Options:
  main arguments: 
    -fa, --fasta PATH             reference sequences in fasta format
                                  [required]
    -fq1, --fastq1 PATH           fastq sequencing file of mate 1  [required]
    -fq2, --fastq2 PATH           fastq sequencing file of mate 2
  common options: 
    --dot_bracket PATH            A csv formatted file that contains dot
                                  bracket info for each sequence
    -pf, --param-file PATH        A yml formatted file to specify parameters
    -ow, --overwrite              overwrites previous results, if not set will
                                  keep previous calculation checkpoints
    -ll, --log-level TEXT         set log level (INFO|WARN|DEBUG|ERROR|FATAL)
    -rob, --restore_org_behavior  retores the original behavior of dreem upon
                                  first release
  map options: 
    --map-overwrite               overwrite mapping calculation
    --skip                        do not perform sequence mapping, not
                                  recommended
    --skip_fastqc                 do not run fastqc for quality control of
                                  sequence data
    --skip_trim_galore            do not run trim galore to remove adapter
                                  sequences at ends
    --tg_q_cutoff TEXT            TODO
    --bt2_alignment_args TEXT     TODO
  bv options: 
    --skip                        skip bit vector generation step, not
                                  recommended
    --bv-overwrite                overwrite bit vector calculation
    --qscore_cutoff TEXT          quality score of read nucleotide, sets to
                                  ambigious if under this val
    --num_of_surbases TEXT        number of bases around a mutation
    --map_score_cutoff TEXT       map alignment score cutoff for a read, read
                                  is discarded if under this value
    --mutation_count_cutoff TEXT  maximum number of mutations in a read
                                  allowable
    --percent_length_cutoff TEXT  read is discarded if less than this percent
                                  of a ref sequence is included
  --help                          Show this message and exit.
```

### running via docker 

If running through docker you can use all the same arguments but use the exe `dreem-docker` 

## Test case

There is example data within the `dreem/test` directory

```shell
cd dreem/
dreem -fa test/resources/case_1/test.fasta -fq1 test/resources/case_1/test_mate1.fastq -fq2 test/resources/case_1/test_mate2.fastq
```

successful output looks as follows

```shell
[16:08 bit_vector.py main] INFO ran at commandline as: 
[16:08 bit_vector.py main] INFO /Users/josephyesselman/installs/anaconda3/envs/py3/bin/dreem -fa test/resources/case_1/test.fasta -fq1 test/resources/case_1/test_mate1.fastq -fq2 test/resources/case_1/test_mate2.fastq
[16:08 bit_vector.py validate_files] INFO fasta file: test/resources/case_1/test.fasta exists
[16:08 bit_vector.py validate_files] INFO fastq file: test/resources/case_1/test_mate1.fastq exists
[16:08 bit_vector.py validate_files] INFO fastq2 file: test/resources/case_1/test_mate2.fastq exists
[16:08 bit_vector.py validate_files] INFO two fastq files supplied, thus assuming paired reads
[16:08 bit_vector.py build_directories] INFO building directory structure
[16:08 mapper.py __init__] INFO bowtie2 2.4.1 detected!
[16:08 mapper.py __init__] INFO fastqc v0.11.9 detected!
[16:08 mapper.py __init__] INFO trim_galore 0.6.6 detected!
[16:08 mapper.py __init__] INFO cutapt 1.18 detected!
[16:08 mapper.py __run_command] INFO running fastqc
[16:08 mapper.py __run_command] INFO fastqc ran without errors
[16:08 mapper.py __run_command] INFO running trim_galore
[16:08 mapper.py __run_command] INFO trim_galore ran without errors
[16:08 mapper.py __run_command] INFO running bowtie2-build
[16:08 mapper.py __run_command] INFO bowtie2-build ran without errors
[16:08 mapper.py __run_command] INFO running bowtie2 alignment
[16:08 mapper.py __run_command] INFO bowtie2 alignment ran without errors
[16:08 mapper.py __run_bowtie_alignment] INFO results for bowtie alignment: 
2500 reads; of these:
  2500 (100.00%) were paired; of these:
    168 (6.72%) aligned concordantly 0 times
    2331 (93.24%) aligned concordantly exactly 1 time
    1 (0.04%) aligned concordantly >1 times
93.28% overall alignment rate
[16:08 mapper.py __run_picard_bam_convert] INFO Converting BAM file to SAM file format
[16:08 mapper.py __run_command] INFO running picard BAM conversion
[16:08 mapper.py __run_command] INFO picard BAM conversion ran without errors
[16:08 mapper.py __run_picard_sort] INFO sorting BAM file
[16:08 mapper.py __run_command] INFO running picard BAM sort
[16:08 mapper.py __run_command] INFO picard BAM sort ran without errors
[16:08 mapper.py run] INFO finished mapping!
[16:08 sam.py run] INFO starting bitvector generation
[16:08 sam.py __run_command] INFO running picard SAM convert
[16:08 sam.py __run_command] INFO picard SAM convert ran without errors
[16:08 sam.py run] INFO MUTATION SUMMARY
| name          |   reads |   aligned |   no_mut |   1_mut |   2_mut |   3_mut |   3plus_mut |
|---------------|---------|-----------|----------|---------|---------|---------|-------------|
| mttr-6-alt-h3 |    2332 |      2331 |    46.42 |   36.81 |   13.21 |    3.05 |        0.04 |

```

The same command can be run through docker

```shell
cd dreem/
dreem-docker -fa test/resources/case_1/test.fasta -fq1 test/resources/case_1/test_mate1.fastq -fq2 test/resources/case_1/test_mate2.fastq
```

output as follows 

```shell
[16:12 RUN_DOCKER main] INFO DOCKER CMD:
docker run -v /Users/josephyesselman/projects/dreem/test/resources/case_1/test.fasta:/data/test.fasta -v /Users/josephyesselman/projects/dreem/test/resources/case_1/test_mate1.fastq:/data/test_mate1.fastq -v /Users/josephyesselman/projects/dreem/test/resources/case_1/test_mate2.fastq:/data/test_mate2.fastq --name dreem_cont -t dreem  dreem -fa test.fasta -fq1 test_mate1.fastq -fq2 test_mate2.fastq --log-level INFO 
[21:12 bit_vector.py main] INFO ran at commandline as: 
[21:12 bit_vector.py main] INFO /opt/conda/envs/py3/bin/dreem -fa test.fasta -fq1 test_mate1.fastq -fq2 test_mate2.fastq --log-level INFO
[21:12 bit_vector.py validate_files] INFO fasta file: test.fasta exists
[21:12 bit_vector.py validate_files] INFO fastq file: test_mate1.fastq exists
[21:12 bit_vector.py validate_files] INFO fastq2 file: test_mate2.fastq exists
[21:12 bit_vector.py validate_files] INFO two fastq files supplied, thus assuming paired reads
[21:12 bit_vector.py build_directories] INFO building directory structure
[21:12 mapper.py __init__] INFO bowtie2 2.2.9 detected!
[21:12 mapper.py __init__] INFO fastqc v0.11.4 detected!
[21:12 mapper.py __init__] INFO trim_galore 0.6.6 detected!
[21:12 mapper.py __init__] INFO cutapt 1.18 detected!
[21:12 mapper.py __run_command] INFO running fastqc
[21:12 mapper.py __run_command] INFO fastqc ran without errors
[21:12 mapper.py __run_command] INFO running trim_galore
[21:12 mapper.py __run_command] INFO trim_galore ran without errors
[21:12 mapper.py __run_command] INFO running bowtie2-build
[21:12 mapper.py __run_command] INFO bowtie2-build ran without errors
[21:12 mapper.py __run_command] INFO running bowtie2 alignment
[21:12 mapper.py __run_command] INFO bowtie2 alignment ran without errors
[21:12 mapper.py __run_bowtie_alignment] INFO results for bowtie alignment: 
2500 reads; of these:
  2500 (100.00%) were paired; of these:
    38 (1.52%) aligned concordantly 0 times
    2333 (93.32%) aligned concordantly exactly 1 time
    129 (5.16%) aligned concordantly >1 times
98.48% overall alignment rate
[21:12 mapper.py __run_picard_bam_convert] INFO Converting BAM file to SAM file format
[21:12 mapper.py __run_command] INFO running picard BAM conversion
[21:12 mapper.py __run_command] INFO picard BAM conversion ran without errors
[21:12 mapper.py __run_picard_sort] INFO sorting BAM file
[21:12 mapper.py __run_command] INFO running picard BAM sort
[21:12 mapper.py __run_command] INFO picard BAM sort ran without errors
[21:12 mapper.py run] INFO finished mapping!
[21:12 sam.py run] INFO starting bitvector generation
[21:12 sam.py __run_command] INFO running picard SAM convert
[21:12 sam.py __run_command] INFO picard SAM convert ran without errors
[21:12 sam.py run] INFO MUTATION SUMMARY
| name          |   reads |   aligned |   no_mut |   1_mut |   2_mut |   3_mut |   3plus_mut |
|---------------|---------|-----------|----------|---------|---------|---------|-------------|
| mttr-6-alt-h3 |    2462 |      2333 |    46.55 |   36.69 |    13.2 |    3.04 |        0.04 |
[16:12 RUN_DOCKER main] INFO clean up and copy files from docker
[16:12 RUN_DOCKER main] INFO docker cp dreem_cont:/data/output output
[16:12 RUN_DOCKER main] INFO docker cp dreem_cont:/data/log log
[16:12 RUN_DOCKER main] INFO docker cp dreem_cont:/data/dreem.log dreem.log
[16:12 RUN_DOCKER main] INFO docker rm dreem_cont
dreem_cont

```




