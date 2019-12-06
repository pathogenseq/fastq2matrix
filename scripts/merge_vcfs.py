#! /usr/bin/env python

import sys
import argparse
import subprocess
from fastq2matrix import run_cmd, nofile, nofolder, vcf_class


def main(args):
    FAILED_SAMPLES = open("%s.failed_samples.log" % args.prefix, "w")
    params = {"threads": args.threads, "prefix": args.prefix, "ref": args.ref, "map_file": f"{args.prefix}.map"}

    with open(params["map_file"],"w") as O:
        # Set up list to hold sample names
        samples = []
        # Loop through sample-file and do (1) append samples to list, (2) write sample to map file and (3) check for VCF index
        for line in open(args.sample_file):
            sample = line.rstrip()
            vcf_file = f"{args.vcf_dir}/{sample}{args.vcf_extension}"
            sys.stderr.write(f"Looking for {vcf_file}\n")
            if os.path.isfile(vcf_file):
                sys.stderr.write("...OK\n")
            else:
                sys.stderr.write("...Not found...skipping\n")
                continue
            # filecheck(vcf_file)
            if args.ignore_missing and nofile(vcf_file):
                FAILED_SAMPLES.write("%s\tno_file\n" % sample)
                continue
            if nofile(f"{vcf_file}.validated"):
                if nofile(f"{vcf_file}.tbi"):
                    run_cmd(f"tabix {vcf_file}")
                run_cmd(f"gatk ValidateVariants -R {args.ref} -V {vcf_file} -gvcf && touch {vcf_file}.validated")
                if nofile(f"{vcf_file}.validated"):
                    FAILED_SAMPLES.write("%s\tno_validation\n" % sample)
                    continue
            samples.append(sample)
            O.write("%s\t%s\n" % (sample,vcf_file))
            if nofile(f"{vcf_file}.tbi"):
                run_cmd(f"bcftools index --tbi {vcf_file}")
    stages = {"dbimport":1,"genotype":2,"filtering":3,"fasta":4,"matrix":5,"pca":6}
    # Create .dict file (GATK fasta index) has been created for the reference
    if nofile("%s.dict" % args.ref.replace(".fasta","").replace(".fa","")):
        run_cmd("gatk CreateSequenceDictionary -R %(ref)s" % params)
    # Create .fai file (SAMtools fasta index) has been created for the reference
    if nofile("%s.fai" % args.ref.replace(".fasta","").replace(".fa","")):
        run_cmd("samtools faidx %(ref)s" % params)
    if nofolder("%(prefix)s_genomics_db" % params) or stages[args.redo]<=1:
        run_cmd("gatk GenomicsDBImport --genomicsdb-workspace-path %(prefix)s_genomics_db -L Chromosome --sample-name-map %(map_file)s --reader-threads %(threads)s --batch-size 500" % params, verbose=2)
    if nofile("%(prefix)s.raw.vcf.gz" % params) or stages[args.redo]<=2:
        run_cmd("gatk --java-options \"-Xmx40g\" GenotypeGVCFs -R %(ref)s -V gendb://%(prefix)s_genomics_db -O %(prefix)s.raw.vcf.gz" % params, verbose=2)

parser = argparse.ArgumentParser(description='VCF mergin pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--sample-file',help='sample file',required=True)
parser.add_argument('--prefix',help='Prefix for files',required=True)
parser.add_argument('--ref',help='reference file',required=True)
parser.add_argument('--vcf-dir',default="./vcf/", type=str, help='VCF firectory')
parser.add_argument('--vcf-extension',default=".gatk.vcf.gz", type=str, help='VCF extension')
parser.add_argument('--threads',default=4, type=int, help='Number of threads for parallel operations')
parser.add_argument('--ignore-missing', action="store_true", help='If this option is set, missing samples are ignored')
parser.add_argument('--redo',type=str,choices=["dbimport","genotype","filtering","fasta","matrix","pca"])
parser.add_argument('--no-validate',action="store_true",)
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
