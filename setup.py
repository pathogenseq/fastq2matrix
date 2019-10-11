import setuptools


setuptools.setup(

	name="fastq2matrix",
	version="1.6.1",
	packages=["fastq2matrix",],
	license="GPL3",
	long_description="Pathogen profiling tool",
	scripts=[
		'scripts/fastq2vcf.py',
	]
)
