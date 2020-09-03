#! /usr/bin/env python

import sys
import os
import argparse
import subprocess
from fastq2matrix import run_cmd, cmd_out, nofile, nofolder, vcf_class, get_contigs_from_fai, get_random_file, filecheck, foldercheck
import json
from datetime import date

def main_import(args):
    FAILED_SAMPLES = open("%s.failed_samples.log" % args.prefix, "w")
    params = vars(args)
    params["map_file"]= f"{args.prefix}.map"

    with open(params["map_file"],"w") as O:
        # Set up list to hold sample names
        samples = []
        # Loop through sample-file and do (1) append samples to list, (2) write sample to map file and (3) check for VCF index
        for line in open(args.sample_file):
            sample = line.rstrip()
            vcf_file = f"{args.vcf_dir}/{sample}{args.vcf_extension}"
            sys.stderr.write(f"Looking for {vcf_file}")
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
    # Create .dict file (GATK fasta index) has been created for the reference
    if nofile("%s.dict" % args.ref.replace(".fasta","").replace(".fa","")):
        run_cmd("gatk CreateSequenceDictionary -R %(ref)s" % params)
    # Create .fai file (SAMtools fasta index) has been created for the reference
    if nofile("%s.fai" % args.ref.replace(".fasta","").replace(".fa","")):
        run_cmd("samtools faidx %(ref)s" % params)

    window_cmd = "bedtools makewindows -n %(num_genome_chunks)s -g %(ref)s.fai | awk '{print $1\":\"$2+1\"-\"$3\" \"$1\"_\"$2+1\"_\"$3}'" % params
    if nofile("%(prefix)s.dbconf.json" % params):
        import_cmd = "gatk GenomicsDBImport --genomicsdb-workspace-path %(prefix)s_{2}_genomics_db -L {1} --sample-name-map %(map_file)s --reader-threads %(threads)s --batch-size 500" % params
        run_cmd(f"{window_cmd} | parallel --bar -j {args.threads} --col-sep \" \" {import_cmd}", verbose=2)
        json.dump({"num_genome_chunks":args.num_genome_chunks},open("%(prefix)s.dbconf.json" % params,"w"))
    else:
        conf = json.load(open(filecheck(f"{args.prefix}.dbconf.json")))
        for l in cmd_out(window_cmd):
            row = l.strip().split()
            dirname = "%s_%s_genomics_db" % (args.prefix,row[1])
            sys.stderr.write("Looking for direcotry named %s..." % dirname)
            foldercheck(dirname)
            sys.stderr.write("OK\n")
        import_cmd = "gatk GenomicsDBImport --genomicsdb-update-workspace-path %(prefix)s_{2}_genomics_db -L {1} --sample-name-map %(map_file)s --reader-threads %(threads)s --batch-size 500" % params
        run_cmd(f"{window_cmd} | parallel --bar -j {args.threads} --col-sep \" \" {import_cmd}", verbose=2)

def main_genotype(args):
    conf = json.load(open(filecheck(f"{args.prefix}.dbconf.json")))
    params = vars(args)
    params["num_genome_chunks"] = conf["num_genome_chunks"]
    window_cmd = "bedtools makewindows -n %(num_genome_chunks)s -g %(ref)s.fai | awk '{print $1\":\"$2+1\"-\"$3\" \"$1\"_\"$2+1\"_\"$3}'" % params
    params["window_cmd"] = window_cmd
    # Check folders exist
    for l in cmd_out(window_cmd):
        row = l.strip().split()
        dirname = "%s_%s_genomics_db" % (args.prefix,row[1])
        sys.stderr.write("Looking for direcotry named %s..." % dirname)
        foldercheck(dirname)
        sys.stderr.write("OK\n")

    genotype_cmd = "gatk --java-options \"-Xmx40g\" GenotypeGVCFs -R %(ref)s -V gendb://%(prefix)s_{2}_genomics_db -O %(prefix)s.{2}.genotyped.vcf.gz" % params
    run_cmd(f"{window_cmd} | parallel --bar -j {args.threads} --col-sep \" \" {genotype_cmd}",verbose=2)
    run_cmd("bcftools concat -Oz -o %(prefix)s.%(subfix_vcf)s.genotyped.vcf.gz `%(window_cmd)s | awk '{print \"%(prefix)s.\"$2\".genotyped.vcf.gz\"}'`" % params)
    run_cmd("rm `%(window_cmd)s | awk '{print \"%(prefix)s.\"$2\".genotyped.vcf.gz*\"}'`" % params)


def main_all(args):
    main_import(args)
    main_genotype(args)


parser = argparse.ArgumentParser(description='VCF mergin pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = parser.add_subparsers(help="Task to perform")

#TODO adjust help notes
parser_sub = subparsers.add_parser('all', help='Trim reads using trimmomatic', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--sample-file',help='sample file',required=True)
parser_sub.add_argument('--prefix',help='Prefix for database name',required=True)
parser_sub.add_argument('--ref',help='reference file',required=True)
parser_sub.add_argument('--vcf-dir',default="./vcf/", type=str, help='VCF firectory')
parser_sub.add_argument('--vcf-extension',default=".g.vcf.gz", type=str, help='VCF extension')
parser_sub.add_argument('--threads',default=4, type=int, help='Number of threads for parallel operations')
parser_sub.add_argument('--num-genome-chunks',default=20, type=int, help='Number of chunks to divide the genome into')
parser_sub.add_argument('--ignore-missing', action="store_true", help='If this option is set, missing samples are ignored')
parser_sub.add_argument('--no-validate',action="store_true")
parser_sub.add_argument('--subfix-vcf', default=date.today().strftime('%Y_%m_%d'), type=str, help='Subfix for genotyped vcf')
parser_sub.set_defaults(func=main_all)


parser_sub = subparsers.add_parser('import', help='Trim reads using trimmomatic', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--sample-file',help='sample file',required=True)
parser_sub.add_argument('--prefix',help='Prefix for database name',required=True)
parser_sub.add_argument('--ref',help='reference file',required=True)
parser_sub.add_argument('--vcf-dir',default="./vcf/", type=str, help='VCF firectory')
parser_sub.add_argument('--vcf-extension',default=".g.vcf.gz", type=str, help='VCF extension')
parser_sub.add_argument('--threads',default=4, type=int, help='Number of threads for parallel operations')
parser_sub.add_argument('--num-genome-chunks',default=20, type=int, help='Number of chunks to divide the genome into')
parser_sub.add_argument('--ignore-missing', action="store_true", help='If this option is set, missing samples are ignored')
parser_sub.add_argument('--redo',type=str,choices=["dbimport","genotype","filtering","fasta","matrix","pca"])
parser_sub.add_argument('--no-validate',action="store_true",)
parser_sub.set_defaults(func=main_import)

parser_sub = subparsers.add_parser('genotype', help='Trim reads using trimmomatic', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--prefix',help='Prefix for database name',required=True)
parser_sub.add_argument('--subfix-vcf', default=date.today().strftime('%Y_%m_%d'), type=str, help='Subfix for genotyped vcf')
parser_sub.add_argument('--ref',help='reference file',required=True)
parser_sub.add_argument('--threads',default=4, type=int, help='Number of threads for parallel operations')
parser_sub.set_defaults(func=main_genotype)

args = parser.parse_args()
if vars(args) == {}:
    parser.print_help(sys.stderr)
else:
    args.func(args)
