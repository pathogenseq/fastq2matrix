#! /usr/bin/env python

import sys
import argparse
import fastq2matrix as fm


def main_trim(args):
	fm.run_cmd("trimmomatic PE %(read1)s %(read2)s -baseout %(prefix)s LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36" % vars(args))

def main_map(args):
	fm.run_cmd("bwa mem -t %(threads)s -R \"@RG\\tID:%(prefix)s\\tSM:%(prefix)s\\tPL:Illumina\" %(ref)s %(prefix)s_1P %(prefix)s_2P | samtools view -@ %(threads)s -b - | samtools sort -@ %(threads)s -o %(prefix)s.bam -" % vars(args))
	fm.run_cmd("samtools index -@ %(threads)s %(prefix)s.bam" % vars(args))
	fm.run_cmd("rm %(prefix)s_1P %(prefix)s_2P %(prefix)s_1U %(prefix)s_2U" % vars(args))
	if args.bqsr_vcf:
		fm.run_cmd("gatk BaseRecalibrator -R %(ref)s -I %(prefix)s.bam --known-sites %(bqsr_vcf)s -O %(prefix)s.recal_data.table" % vars(args))
		fm.run_cmd("gatk ApplyBQSR -R %(ref)s -I %(prefix)s.bam --bqsr-recal-file %(prefix)s.recal_data.table -O %(prefix)s.bqsr.bam" % vars(args))
		fm.run_cmd("samtools index -@ %(threads)s %(prefix)s.bqsr.bam" % vars(args))
		fm.run_cmd("rm %(prefix)s.bam" % vars(args))

def main_gatk(args):
	if not args.prefix:
		args.prefix = args.bam.replace(".bam","")
	fm.run_cmd("gatk HaplotypeCaller -I %(bam)s -R %(ref)s -O %(prefix)s.g.vcf.gz -ERC BP_RESOLUTION" % vars(args))
	fm.run_cmd("gatk ValidateVariants -V %(prefix)s.g.vcf.gz -gvcf -R %(ref)s && touch %(prefix)s.g.vcf.gz.validated" % vars(args))

def main_all(args):
	main_trim(args)
	args.bam = args.prefix+".bqsr.bam" if args.bqsr_vcf else args.prefix+".bam"
	main_map(args)
	main_gatk(args)

parser = argparse.ArgumentParser(description='fastq2matrix pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = parser.add_subparsers(help="Task to perform")

parser_sub = subparsers.add_parser('all', help='Trim reads using trimmomatic', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--read1','-1',help='First read file',required=True)
parser_sub.add_argument('--read2','-2',help='Second read file',required=True)
parser_sub.add_argument('--prefix','-p',help='Sample prefix for all results generated',required=True)
parser_sub.add_argument('--ref','-r',help='Second read file',required=True)
parser_sub.add_argument('--threads','-t',default=4,help='Number of threads')
parser_sub.add_argument('--bqsr-vcf','-q',help='VCF file used for bqsr')
parser_sub.set_defaults(func=main_all)

parser_sub = subparsers.add_parser('trim', help='Trim reads using trimmomatic', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--prefix','-p',help='Sample prefix for all results generated',required=True)
parser_sub.set_defaults(func=main_trim)

parser_sub = subparsers.add_parser('map', help='Trim reads using trimmomatic', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--read1','-1',help='First read file',required=True)
parser_sub.add_argument('--read2','-2',help='Second read file',required=True)
parser_sub.add_argument('--prefix','-p',help='Sample prefix for all results generated',required=True)
parser_sub.add_argument('--ref','-r',help='Second read file',required=True)
parser_sub.add_argument('--threads','-t',default=4,help='Number of threads')
parser_sub.add_argument('--bqsr-vcf','-q',help='VCF file used for bqsr')
parser_sub.set_defaults(func=main_map)

parser_sub = subparsers.add_parser('gatk', help='Trim reads using trimmomatic', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--bam','-b',help='First read file',required=True)
parser_sub.add_argument('--ref','-r',help='Second read file',required=True)
parser_sub.add_argument('--prefix','-p',help='Sample prefix for all results generated')
parser_sub.add_argument('--threads','-t',default=4,help='Number of threads')

parser_sub.set_defaults(func=main_gatk)




args = parser.parse_args()
if vars(args) == {}:
	parser.print_help(sys.stderr)
else:
	args.func(args)
