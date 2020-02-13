#! /usr/bin/env python

import sys
import argparse
import fastq2matrix as fm
import os.path


def get_step_num(prefix):
    files = {f"{prefix}.mkdup.bam.bai":1,f"{prefix}.bqsr.bam.bai":2,f"{prefix}.g.vcf.gz.validated":3}
    step = 0
    for f in files:
        if os.path.isfile(f):
            step = files[f]
            sys.stderr.write(f"Found {f}\n")
    return step


def main_trim(args):
    fm.run_cmd("trimmomatic PE -phred33 %(read1)s %(read2)s -baseout %(prefix)s LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36" % vars(args))


def main_map(args):
    args.step = get_step_num(args.prefix)
    if not os.path.isfile(args.ref.replace(".fasta",".fasta.amb")):
        fm.run_cmd("bwa index %s" % args.ref)
    if "trimmed" in vars(args):
        args.reads = "%(prefix)s_1P %(prefix)s_2P" % vars(args)
    else:
        args.reads = "%(read1)s %(read2)s" % vars(args)
    if args.redo or args.step<1:
        fm.run_cmd("bwa mem -t %(threads)s -R \"@RG\\tID:%(prefix)s\\tSM:%(prefix)s\\tPL:Illumina\" %(ref)s %(reads)s | samtools view -@ %(threads)s -b - | samtools fixmate -@ %(threads)s -m - - | samtools sort -@ %(threads)s - | samtools markdup -@ %(threads)s - %(prefix)s.mkdup.bam -" % vars(args))
        if "trimmed" in vars(args):
            fm.run_cmd("rm %(prefix)s_1P %(prefix)s_2P %(prefix)s_1U %(prefix)s_2U" % vars(args))
        fm.run_cmd("samtools index -@ %(threads)s %(prefix)s.mkdup.bam" % vars(args))
        fm.run_cmd("samtools flagstat -@ %(threads)s %(prefix)s.mkdup.bam > %(prefix)s.mkdup.bamstats" % vars(args))
    if args.bqsr_vcf and (args.redo or args.step<2):
        if not os.path.isfile(args.ref.replace(".fasta",".fasta.fai")):
            fm.run_cmd("samtools faidx %s" % args.ref)
        if not os.path.isfile(args.ref.replace(".fasta",".dict")):
            fm.run_cmd("gatk CreateSequenceDictionary -R %s" % args.ref)
        for s in args.bqsr_vcf.split(","):
            if not os.path.isfile(s + ".tbi"):
                fm.run_cmd("bcftools index -t %s" % s)
        args.bqsr_vcf = " ".join(["--known-sites %s" % s for s in args.bqsr_vcf.split(",")])
        fm.run_cmd("gatk BaseRecalibrator -R %(ref)s -I %(prefix)s.mkdup.bam %(bqsr_vcf)s -O %(prefix)s.recal_data.table" % vars(args))
        fm.run_cmd("gatk ApplyBQSR -R %(ref)s -I %(prefix)s.mkdup.bam --bqsr-recal-file %(prefix)s.recal_data.table -O %(prefix)s.bqsr.bam" % vars(args))
        fm.run_cmd("samtools index -@ %(threads)s %(prefix)s.bqsr.bam" % vars(args))
        fm.run_cmd("samtools flagstat -@ %(threads)s %(prefix)s.bqsr.bam > %(prefix)s.bqsr.bamstats" % vars(args))
        fm.run_cmd("rm %(prefix)s.mkdup.bam*" % vars(args))


def main_gatk(args):
    if not args.prefix:
        args.prefix = args.bam.replace(".bam","")
    fm.run_cmd("gatk HaplotypeCaller -I %(bam)s -R %(ref)s -O %(prefix)s.g.vcf.gz -ERC %(erc)s" % vars(args))
    fm.run_cmd("gatk ValidateVariants -V %(prefix)s.g.vcf.gz -gvcf -R %(ref)s && touch %(prefix)s.g.vcf.gz.validated" % vars(args))


def main_all(args):
    args.step = get_step_num(args.prefix)
    args.bam = args.prefix+".bqsr.bam" if args.bqsr_vcf else args.prefix+".mkdup.bam"
    sys.stderr.write(f"Starting at step {args.step+1}")
    if args.redo or args.step<1:
        main_trim(args)
        args.trimmed = True
    if args.redo or args.step<2:
        main_map(args)
    if args.redo or args.step<3:
        sys.stderr.write("Using %(bam)s as the bam file" % vars(args))
        main_gatk(args)

parser = argparse.ArgumentParser(description='fastq2matrix pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = parser.add_subparsers(help="Task to perform")

parser_sub = subparsers.add_parser('all', help='Trim reads using trimmomatic', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--read1','-1',help='First read file',required=True)
parser_sub.add_argument('--read2','-2',help='Second read file',required=True)
parser_sub.add_argument('--prefix','-p',help='Sample prefix for all results generated',required=True)
parser_sub.add_argument('--ref','-r',help='Second read file',required=True)
parser_sub.add_argument('--threads','-t',default=4,help='Number of threads'
