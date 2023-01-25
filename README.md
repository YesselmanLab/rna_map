


# RNA MAP


## How to install


## How to use 

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

