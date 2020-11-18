import setuptools
import os

setuptools.setup(

    name="fastq2matrix",
    version="0.1.3",
    packages=["fastq2matrix"],
    license="MIT",
    long_description="Utilities to get from fastq files to a variant matrix",
    scripts= ["scripts/%s" % x for x in os.listdir("scripts")],

)
