import fastq2matrix as fm
import argparse
import sys

def main(args):
    if args.prefix:
        individual_bams = ["%s/%s%s" % (args.dir,run,args.suffix) for run in args.prefix.split("_")]
        new_id = args.new_id if args.new_id else args.prefix
    elif args.bams:
        individual_bams = args.bams.split(",")
        new_id = args.new_id if args.new_id else "_".join([bam.split("/")[-1].replace(args.suffix,"") for bam in individual_bams])
    elif (not args.prefix and not args.bams) or (args.prefix and args.bams):
        sys.stderr.write("Need wither '--bams' or '--prefix'... Exiting!\n")
        quit()
    if len(individual_bams)==1:
        sys.stderr.write("Need more than one bam... Exiting!\n")
        quit()
    for bam in individual_bams:
        fm.filecheck(bam)
    new_bamfile = "%s/%s%s" % (args.dir,new_id,args.suffix)
    tmp_bamfile = fm.get_random_file()
    tmp_file = fm.get_random_file()
    with open(tmp_file,"w") as O:
        for l in fm.cmd_out("samtools view -H %s" % individual_bams[0]):
            row = l.strip().split("\t")
            if row[0]=="@RG":
                continue
                row[1] = "ID:%s" % new_id
                row[2] = "SM:%s" % new_id
            O.write("%s\n" % "\t".join(row))

    fm.run_cmd("samtools merge -@ %s - %s | samtools reheader -i %s - | samtools addreplacerg -@ %s - -r 'ID:%s\\tSM:%s\\tPL:Illumina' -o %s" % (
        args.threads," ".join(individual_bams), tmp_file, args.threads,new_id, new_id, new_bamfile)
    )
    fm.run_cmd("samtools index %s" % new_bamfile)
    fm.rm_files([tmp_file,tmp_bamfile])

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--bams',help='Comma seperated BAM file list')
parser.add_argument('--new-id',help='New id for BAM')
parser.add_argument('--prefix',help='Prefix')
parser.add_argument('--suffix',default=".bqsr.bam",help='BAM file suffix')
parser.add_argument('--dir',default=".",help='Directory containing bams')
parser.add_argument('--threads',default=4, type=int, help='Number of threads for parallel operations')
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
