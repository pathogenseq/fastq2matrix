import sys
import argparse
import fastq2matrix as fm

def main(args):
    args.region_arg=""

    if args.variant_caller=="gatk":
        if args.bed:
            args.region_arg = "-L %s" % args.bed
        fm.run_cmd("gatk HaplotypeCaller -R %(ref)s %(region_arg)s -I %(bam)s -O %(out)s.vcf.gz" % vars(args))
    elif args.variant_caller=="bcftools":
        if args.bed:
            args.region_arg = "-R %s" % args.bed
        fm.run_cmd("bcftools mpileup -f %(ref)s %(region_arg)s %(bam)s | bcftools call -mv -Oz -o %(out)s.vcf.gz" % vars(args))
    elif args.variant_caller=="freebayes":
        if args.bed:
            args.region_arg = "-t %s" % args.bed
        fm.run_cmd("freebayes -f %(ref)s %(region_arg)s %(bam)s | bgzip -c > %(out)s.vcf.gz" % vars(args))
    else:
        quit("Unknown variant caller! Exiting!")

    fm.run_cmd("tabix -f %(out)s.vcf.gz" % vars(args))

    if args.bed:
        fm.run_cmd("bedtools coverage -a %(bed)s -b %(bam)s -d | awk '$NF<%(depth_cutoff)s {print $1\"\\t\"$2+$(NF-1)-2\"\\t\"$2+$(NF-1)-1}' > %(out)s.depth_mask.bed" % vars(args))
    else:
        fm.run_cmd("bedtools genomecov -ibam %(bam)s  -d | awk '$NF<%(depth_cutoff)s {print $1\"\\t\"$2-1\"\\t\"$2}' > %(out)s.depth_mask.bed" % vars(args))


    for l in fm.cmd_out("wc -l %(out)s.depth_mask.bed" % vars(args)):
        num_lines = int(l.strip().split()[0])

    args.mask_arg = "-m %(out)s.depth_mask.bed -M N" % vars(args) if num_lines>0 else ""

    region_names = {}
    if args.bed:
        regions_file = args.out+".regions.txt"
        with open(regions_file, "w") as O:
            for l in open(args.bed):
                row = l.strip().split()
                r = "%s:%s-%s" % (row[0],row[1],row[2])
                O.write(r+"\n")
                if len(row)>3:
                    region_names[r] = row[3]

        args.region_arg = "-r %s" % regions_file
        consensus_cmd = "samtools faidx %(ref)s %(region_arg)s | bcftools consensus %(out)s.vcf.gz %(mask_arg)s" % vars(args)
    else:
        consensus_cmd = "bcftools consensus -f %(ref)s %(out)s.vcf.gz %(mask_arg)s" % vars(args)

    with open(args.out+".consensus.fa","w") as O:
        for l in fm.cmd_out(consensus_cmd):
            if l[0]==">":
                r = l.strip()[1:]
                O.write(">%s %s\n" % (args.out,region_names.get(r,r)))
            else:
                O.write(l+"\n")

parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--bam',type=str,help='Bam file',required=True)
parser.add_argument('--ref',type=str,help='Reference file',required=True)
parser.add_argument('--out',type=str,help='Output prefix',required=True)
parser.add_argument('--bed',type=str,help='Bed file with regions to be extracted (optional 4th column to contain region name which will be present in consensus file)')
parser.add_argument('--variant-caller',choices=["gatk","freebayes","bcftools"],default="freebayes",type=str,help='Variant caller to be used')
parser.add_argument('--depth-cutoff',default=10,type=int,help='Sites with less than this cutoff will be marked as missing in the consensus file')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
