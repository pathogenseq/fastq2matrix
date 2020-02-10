#! /usr/bin/env python

import sys
import os
import argparse
import subprocess
from fastq2matrix import run_cmd, nofile, nofolder, vcf_class, get_contigs_from_fai


def main(args):
    params = {"threads": args.threads, "prefix": args.prefix, "ref": args.ref, "map_file": f"{args.prefix}.map", "merged_file": args.merged_file, "include": args.include_regions, "vqslod": args.vqslod, "miss": args.missing_sample_cutoff, "mix":args.cutoff_mix_GT, "gff_file": args.gff_file}
        
    if args.include_regions:
        if not os.path.isfile("%(merged_file)s.tbi" % params):
            run_cmd("bcftools index -t %(merged_file)s" % params)
        params["vcf_in"] = params["merged_file"].replace(".genotyped.vcf.gz",".in.genotyped.vcf.gz")
        run_cmd("bcftools view -R %(include)s -O z -o %(vcf_in)s %(merged_file)s" % params)
        run_cmd("bcftools index -t %(vcf_in)s" % params)
        params["merged_file"] = params["vcf_in"]
    if not os.path.isfile(args.ref.replace(".fasta",".dict")):
        run_cmd("gatk CreateSequenceDictionary -R %s" % args.ref)
    for s in args.bqsr_vcf.split(","):
        if not os.path.isfile(s + ".tbi"):
            run_cmd("bcftools index -t %s" % s)
    if not os.path.isfile("%(merged_file)s.tbi" % params):
            run_cmd("bcftools index -t %(merged_file)s" % params)
    params["bqsr_vcf_mer"] = " ".join(["--resource:pf_crosses,known=false,training=true,truth=true,prior=15.0 %s " % s for s in args.bqsr_vcf.split(",")])
    params["output"] = params["merged_file"].replace(".genotyped.vcf.gz",".recal")
    ## Calculating calibration model
    run_cmd("gatk VariantRecalibrator -R %(ref)s -V %(merged_file)s %(bqsr_vcf_mer)s -an QD -an FS -an SOR -an DP --max-gaussians 8 --mq-cap-for-logit-jitter-transform 70 -mode SNP -O %(prefix)s.snps.recal --tranches-file %(prefix)s.snps.tranches --rscript-file %(prefix)s.snps.plots.R" % params)
    run_cmd("gatk VariantRecalibrator -R %(ref)s -V %(merged_file)s %(bqsr_vcf_mer)s -an QD -an DP -an SOR -an FS --max-gaussians 4 --mq-cap-for-logit-jitter-transform 70 -mode INDEL -O %(prefix)s.indel.recal --tranches-file %(prefix)s.indel.tranches --rscript-file %(prefix)s.indel.plots.R" % params)
    ## Applying calibration model and obtaining VQSLOD
    run_cmd("gatk ApplyVQSR -R %(ref)s -V %(merged_file)s -O %(prefix)s.vqslod.snps.vcf.gz --truth-sensitivity-filter-level 99.0 --tranches-file %(prefix)s.snps.tranches --recal-file %(prefix)s.snps.recal -mode SNP" % params)
    run_cmd("gatk ApplyVQSR -R %(ref)s -V %(merged_file)s -O %(prefix)s.vqslod.indel.vcf.gz --truth-sensitivity-filter-level 99.0 --tranches-file %(prefix)s.indel.tranches --recal-file %(prefix)s.indel.recal -mode INDEL" % params)
    ## Filtering based on VQSLOD
    run_cmd("bcftools view -i 'VQSLOD>%(vqslod)s' -O z -o %(prefix)s.vqslod.filt.snps.vcf.gz  %(prefix)s.vqslod.snps.vcf.gz" % params)
    run_cmd("bcftools view -i 'VQSLOD>%(vqslod)s' -O z -o %(prefix)s.vqslod.filt.indel.vcf.gz  %(prefix)s.vqslod.indel.vcf.gz" % params)
    ## Annotating filtered files VQSLOD
    run_cmd("bcftools index -t %(prefix)s.vqslod.filt.snps.vcf.gz" % params)
    run_cmd("bcftools index -t %(prefix)s.vqslod.filt.indel.vcf.gz" % params)
    ## Add sample filtering by missing
    if params["miss"] == "0":
        run_cmd("mv %(prefix)s.vqslod.filt.snps.vcf.gz %(prefix)s.miss%(miss)s.vqslod.filt.snps.vcf.gz" % params)
    else:
        run_cmd("plink --vcf %(prefix)s.vqslod.filt.snps.vcf.gz --mind %(miss)s --recode vcf --allow-extra-chr --out %(prefix)s_plink" % params)
        run_cmd("grep -P \"^#CHROM\" %(prefix)s_plink.vcf | awk '{ $1=\"\"; $2=\"\";$3=\"\"; $4=\"\";$5=\"\"; $6=\"\";$7=\"\"; $8=\"\";$9=\"\"; print}' | sed 's/ /\\n/g' | tail -n+10 > %(prefix)s_new" % params)
        run_cmd("bcftools view -S %(prefix)s_new --threads 20 -O z -o  %(prefix)s.miss%(miss)s.vqslod.filt.snps.vcf.gz %(prefix)s.vqslod.filt.snps.vcf.gz" % params)
    ## Add set GT
    run_cmd("bcftools view %(prefix)s.miss%(miss)s.vqslod.filt.snps.vcf.gz | setGT.py --fraction %(mix)s | bcftools view -O z -c 1 -o %(prefix)s.GT.miss%(miss)s.vqslod.filt.snps.vcf.gz" % params)
    ## Select only biallelic
    run_cmd("bcftools view -m2 -M2 %(prefix)s.GT.miss%(miss)s.vqslod.filt.snps.vcf.gz --threads 20 -O z -o %(prefix)s.bi.GT.miss%(miss)s.vqslod.filt.snps.vcf.gz" % params)
    ## Add CSQ annotation
    run_cmd("bcftools csq -p m -f %(ref)s -g %(gff_file)s %(prefix)s.bi.GT.miss%(miss)s.vqslod.filt.snps.vcf.gz -O z -o %(prefix)s.csq.bi.GT.miss%(miss)s.vqslod.filt.snps.vcf.gz" % params)

parser = argparse.ArgumentParser(description='Filtering merged VCF using GATK VQSR',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--merged-file',help='Merged VCF file',required=True)
parser.add_argument('--prefix',help='Prefix for output files',required=True)
parser.add_argument('--ref',help='Reference file',required=True)
parser.add_argument('--bqsr-vcf','-q',help='VCF file used for bqsr, if multiple files available use comma separated list',required=True)
parser.add_argument('--include-regions', help='BED file with regions to be included')
parser.add_argument('--threads',default=4, type=int, help='Number of threads for parallel operations')
parser.add_argument('--vqslod',default=0, type=int, help='VQSLOD filter to exclude SNPs, any SNP wiht VQSLOD< threshold will be removed')
parser.add_argument('--missing-sample-cutoff',default="0", help='Threshold proportion for sample filtering based on SNP positions missing, samples with proportion of SNPs missing higher than threshold will be excluded')
parser.add_argument('--cutoff-mix-GT',default="0.8",help='Threshold for heterozygous GT assignment, SNPs with MAF of second allele higher than threshold will be relabelled as heterozygous calls')
parser.add_argument('--gff-file', help='GFF file for SNP annotation')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
