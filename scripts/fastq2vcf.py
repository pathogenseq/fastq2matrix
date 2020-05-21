#! /usr/bin/env python

import sys
import argparse
import fastq2matrix as fm
import os.path


def get_step_num(prefix):
    files = {f"{prefix}.mkdup.bam.bai":1,f"{prefix}.bqsr.bam.bai":2,f"{prefix}.bqsr.cram.crai":2,f"{prefix}.mkdup.cram.crai":2,f"{prefix}.g.vcf.gz.validated":3}
    step = 0
    for f in files:
        if os.path.isfile(f):
            step = files[f]
            sys.stderr.write(f"Found {f}\n")
    return step


def main_trim(args):
    if args.single:
        fm.run_cmd("trimmomatic SE -phred33 %(read1)s %(prefix)s_trimmed.fq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36" % vars(args))
    else:
        fm.run_cmd("trimmomatic PE -phred33 %(read1)s %(read2)s -baseout %(prefix)s LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36" % vars(args))


def main_map(args):
    args.step = get_step_num(args.prefix)



    if "trimmed" in vars(args) and args.single:
        args.reads = "%(prefix)s_trimmed.fq" % vars(args)
    elif "trimmed" in vars(args) and not args.single:
        args.reads = "%(prefix)s_1P %(prefix)s_2P" % vars(args)
    elif "trimmed" not in vars(args) and args.single:
        args.reads = "%(read1)s %(read2)s" % vars(args)
    elif "trimmed" not in vars(args) and not args.single:
        args.reads = "%(read1)s" % vars(args)
    if args.redo or args.step<1:
        fm.run_cmd("bwa mem -t %(threads)s -R \"@RG\\tID:%(prefix)s\\tSM:%(prefix)s\\tPL:Illumina\" %(ref)s %(reads)s | samtools view -@ %(threads)s -b - | samtools fixmate -@ %(threads)s -m - - | samtools sort -@ %(threads)s - | samtools markdup -@ %(threads)s - %(prefix)s.mkdup.bam -" % vars(args))
        if "trimmed" in vars(args) and args.single:
            fm.run_cmd("rm %(reads)s" % vars(args))
        if "trimmed" in vars(args) and not args.single:
            fm.run_cmd("rm %(prefix)s_1P %(prefix)s_2P %(prefix)s_1U %(prefix)s_2U" % vars(args))
        fm.run_cmd("samtools index -@ %(threads)s %(prefix)s.mkdup.bam" % vars(args))
        fm.run_cmd("samtools flagstat -@ %(threads)s %(prefix)s.mkdup.bam > %(prefix)s.mkdup.bamstats" % vars(args))
    if args.bqsr_vcf and (args.redo or args.step<2):
        for vcf in args.bqsr_vcf.split(","):
            fm.tabix_vcf(vcf)
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

def convert_to_cram(bam_file,ref_file,threads):
    cram_file = bam_file.replace(".bam",".cram")
    fm.run_cmd("samtools view -@ %s -C %s -o %s -T %s" %(threads,bam_file,cram_file,ref_file))
    fm.run_cmd("samtools index %s" % cram_file)
    fm.run_cmd("rm %s %s.bai" % (bam_file,bam_file))

def main_all(args):
    fm.create_seq_dict(args.ref)
    fm.bwa_index(args.ref)
    fm.faidx(args.ref)

    args.step = get_step_num(args.prefix)

    if not args.read2 and not args.single:
        sys.stderr.write("Second read is not provided, please check... Exiting!\n")
        quit()

    args.bam = args.prefix+".bqsr.bam" if args.bqsr_vcf else args.prefix+".mkdup.bam"

    sys.stderr.write("Starting at step %s\n" % (args.step+1))
    if args.redo or args.step<1:
        main_trim(args)
        args.trimmed = True
    if args.redo or args.step<2:
        main_map(args)
    if args.redo or args.step<3:
        sys.stderr.write("Using %(bam)s as the bam file" % vars(args))
        main_gatk(args)
    if args.cram:
        if args.redo or args.step<4:
            convert_to_cram(args.bam,args.ref,args.threads)

parser = argparse.ArgumentParser(description='fastq2matrix pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = parser.add_subparsers(help="Task to perform")

parser_sub = subparsers.add_parser('all', help='Trim reads using trimmomatic', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--read1','-1',help='First read file',required=True)
parser_sub.add_argument('--read2','-2',help='Second read file')
parser_sub.add_argument('--prefix','-p',help='Sample prefix for all results generated',required=True)
parser_sub.add_argument('--ref','-r',help='Second read file',required=True)
parser_sub.add_argument('--threads','-t',default=4,help='Number of threads')
parser_sub.add_argument('--bqsr-vcf','-q',help='VCF file used for bqsr')
parser_sub.add_argument('--erc',default="GVCF", choices=["GVCF","BP_RESOLUTION"], help='Choose ERC type on GATK')
parser_sub.add_argument('--redo',action="store_true",help='Redo everything')
parser_sub.add_argument('--single',action="store_true",help='Redo everything')
parser_sub.add_argument('--cram',action="store_true",help='Redo everything')
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
parser_sub.add_argument('--redo',action="store_true",help='Redo everything')
parser_sub.set_defaults(func=main_map)

parser_sub = subparsers.add_parser('gatk', help='Trim reads using trimmomatic', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--bam','-b',help='First read file',required=True)
parser_sub.add_argument('--ref','-r',help='Second read file',required=True)
parser_sub.add_argument('--prefix','-p',help='Sample prefix for all results generated')
parser_sub.add_argument('--erc',default="GVCF", choices=["GVCF","BP_RESOLUTION"], help='Choose ERC type on GATK')
parser_sub.add_argument('--threads','-t',default=4,help='Number of threads')

parser_sub.set_defaults(func=main_gatk)


args = parser.parse_args()
if vars(args) == {}:
    parser.print_help(sys.stderr)
else:
    args.func(args)
