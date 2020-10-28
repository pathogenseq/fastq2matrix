#! /usr/bin/env python
import sys
import subprocess
import argparse
import os
from fastq2matrix import vcf_class, nofile


def main(args):
	if nofile(args.vcf): quit("Can't find %s... Exiting!" % args.vcf)
	vcf = vcf_class(args.vcf)
	vcf.vcf_to_matrix(args.no_iupacgt,na=args.na)

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--vcf',help='VCF file',required=True)
parser.add_argument('--no-iupacgt', dest='no_iupacgt', action='store_true')
parser.add_argument('--na', default="N",help = "Character to use for missing data value")
parser.add_argument('--threads',default=4, type=int, help='Number of threads for parallel operations')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
