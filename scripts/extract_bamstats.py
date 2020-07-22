import sys
import argparse
import os
from tqdm import tqdm


#15650507 + 0 in total (QC-passed reads + QC-failed reads)
#0 + 0 secondary
#34633 + 0 supplementary
#166695 + 0 duplicates
#15546150 + 0 mapped (99.33% : N/A)
#15615874 + 0 paired in sequencing
#7807937 + 0 read1
#7807937 + 0 read2
#14991678 + 0 properly paired (96.00% : N/A)
#15499656 + 0 with itself and mate mapped
#11861 + 0 singletons (0.08% : N/A)
#0 + 0 with mate mapped to a different chr
#0 + 0 with mate mapped to a different chr (mapQ>=5)

def main(args):
    if args.sample_file:
        samples = [x.rstrip() for x in open(args.sample_file).readlines()]
    else:
        samples = [x.replace(".bqsr.bamstats","") for x in os.listdir("%s/" % args.dir) if x[-14:]==".bqsr.bamstats"]
    for s in tqdm(samples):
        res = []
        for i,l in enumerate(open(f"{args.dir}/{s}.bqsr.bamstats")):
            row = l.rstrip().split()
            if i==4:
                res.append(str(row[0]))
                res.append(str(row[4][1:-1]))
        if args.depth:
            res.append(open(f"{args.dir}/{s}.median_depth").readline().strip())
        print("%s\t%s" % (s,"\t".join(res)))

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--sample-file',type=str,help='Sample file')
parser.add_argument('--dir',default=".",type=str,help='Directory')
parser.add_argument('--depth',action="store_true",help='Add depth info')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
