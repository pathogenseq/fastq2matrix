#! /usr/bin/env python
import sys
from tqdm import tqdm
import argparse
import re
def main(args):
    ad_cutoff = args.fraction
    for l in tqdm(sys.stdin):
        something_changed = False
        row = l.strip().replace("|","/").split()
        if l[0]=="#":
            sys.stdout.write(l)
            continue
        alleles = [row[3]]+row[4].split(",")
        if len(alleles)>9: continue

        if "*" in row[4]:
            something_changed = True
            for i in range(9,len(row)):
                if row[i][0]=="." or row[i][2]==".": continue
                if alleles[int(row[i][0])]=="*":
                    row[i][0] = "."
                if alleles[int(row[i][2])]=="*":
                    row[i][2] = "."
        uniq_mixed_genotypes = set([x for x in re.findall("[0-9]/[0-9]",l) if x[0]!=x[2]])
        if len(uniq_mixed_genotypes)>1:
            for i in range(9,len(row)):
                if row[i][:3] not in uniq_mixed_genotypes: continue
                fmt = row[i].split(":")
                ad = [int(x) for x in fmt[1].split(",")]
                total_ad = sum(ad)
                if total_ad==0:continue
                adf = [ad[j]/total_ad for j in range(len(ad))]
                if max(adf)>=ad_cutoff:
                    new_gt = adf.index(max(adf))
                    fmt[0] = f"{new_gt}/{new_gt}"
                    something_changed = True
                    row[i] = ":".join(fmt)
        if something_changed:
            sys.stdout.write("\t".join(row)+"\n")
        else:
            sys.stdout.write(l)

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--fraction',default=0.7,type=float,help='Fraction of coverage to assign major')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
