import setuptools


setuptools.setup(

    name="fastq2matrix",
    version="0.1.3",
    packages=["fastq2matrix"],
    license="MIT",
    long_description="Utilities to get from fastq files to a variant matrix",
    scripts= [
        'scripts/fastq2vcf.py',
        'scripts/merge_vcfs.py',
        'scripts/setGT.py',
        'scripts/vcf2fasta.py',
		'scripts/vcf2dist.py',
        'scripts/filter_merged_vcf.py',
        'scripts/vcf2matrix.py',
        'scripts/merge_bams.py',
        'scripts/ancestral_reconstruction.py',
        'scripts/filter_tb_vcf.py',
        'scripts/extract_bamstats.py'
        ],

)
