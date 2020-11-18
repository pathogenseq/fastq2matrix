# fastq2matrix

This repo contains scripts to help go from fastq files to merged variant datatype (vcf, fasta, matrix).
The name comes from the old days when we used to rely on a text based variant matrix file.

## Installation

```
conda install python=3.7 bwa samtools bcftools parallel datamash gatk4=4.1.4.1 delly tqdm trimmomatic minimap2 biopython bedtools r-ggplot2 iqtree

git clone https://github.com/pathogenseq/fastq2matrix.git
python setup.py install
```

## Updating your install
If you are not up to date with the latest version of this package then you can run the following code to update. Just remember to replace `/path/to/fastq2matrix` with the location where you initially cloned the repo.
```
cd /path/to/fastq2matrix
git pull
python setup.py install
```

## More info
For an example of how to go from fastq files to a multi-sample vcf file check out [this page](https://jodyphelan.gitbook.io/tutorials/ngs/fastq-to-vcf).
