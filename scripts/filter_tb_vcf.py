#! /usr/bin/env python
import sys
import argparse
import fastq2matrix as fm


def main(args):
    original_vcf = fm.vcf_class(args.vcf,threads=args.threads)
    args.prefix = original_vcf.prefix
    args.filename = original_vcf.filename
    args.indels_cmd = "" if args.keep_indels else "bcftools view -V indels | "
    args.exclude_cmd = f"bcftools view -T ^{args.exclude_bed} |"  if args.exclude_bed else ""
    args.annotation_drop = "" if args.keep_genotype_info else "bcftools annotate -x ^FORMAT/GT | "

    args.window_cmd = "bedtools makewindows -g %(ref)s.fai -n 20 | awk '{print $1\":\"$2+1\"-\"$3\" \"$1\"_\"$2+1\"_\"$3}'" % vars(args)
    args.filter_cmd = (
        "%(indels_cmd)s"
        "%(exclude_cmd)s"
        "setGT.py | "
        "%(annotation_drop)s"
        "bcftools view -c 1 -a -Ou | "
        "bcftools filter -e 'GT=\\\"het\\\"' -S . | "
        "bcftools view -i 'F_PASS(GT!=\\\"mis\\\")>%(site_missing)s' | "
        "bcftools view -c 1 | "
        "bcftools +fill-tags | "
        "bcftools view -e 'AF==1 || AF==0' | "
        "bcftools norm -f %(ref)s" % vars(args)
    )
    if args.keep_indels:
        args.final_file = "%(prefix)s.filtered.vcf.gz" % vars(args)
    else:
        args.final_file = "%(prefix)s.filtered_no_indels.vcf.gz" % vars(args)

    fm.run_cmd("%(window_cmd)s | parallel -j %(threads)s --col-sep \" \" \"bcftools view  %(filename)s -r {1} | %(filter_cmd)s > %(prefix)s.{2}.tmp.txt\"" % vars(args))
    fm.run_cmd("bcftools concat -Oz -o %(final_file)s `%(window_cmd)s | awk '{print \"%(prefix)s.\"$2\".tmp.txt\"}'`" % vars(args))
    fm.run_cmd("rm `%(window_cmd)s | awk '{print \"%(prefix)s.\"$2\".tmp.txt*\"}'`" % vars(args))


parser = argparse.ArgumentParser(description='XXX pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--vcf',help='VCF file',required=True)
parser.add_argument('--ref',help='Reference fasta file',required=True)
parser.add_argument('--site-missing',default=0.9,type=float,help='Minimum non fraction of non-missing data to keep site')
parser.add_argument('--threads',default=2,type=int,help='Number of threads to use')
parser.add_argument('--keep-indels',action="store_true",help='Keep in output')
parser.add_argument('--exclude-bed',help='BED file with regions to remove')
parser.add_argument('--keep-genotype-info',action="store_true",help='Keep information like depth and genotype probabilities for each call')
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
