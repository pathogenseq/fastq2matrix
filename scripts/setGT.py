#! /usr/bin/env python
import sys
from tqdm import tqdm
import argparse

def main(args):
    ad_cutoff = args.fraction
    for l in tqdm(sys.stdin):
    	something_changed = False
    	row = l.strip().split()
    	if l[0]=="#":
    		sys.stdout.write(l)
    		continue
    	for i in range(9,len(row)):
    		fmt = row[i].split(":")
    		ad = [int(x) for x in fmt[1].split(",")]
    		total_ad = sum(ad)
    		if total_ad==0:continue
    		adf = [ad[j]/total_ad for j in range(len(ad))]
    		if max(adf)>=ad_cutoff:
    			gt = adf.index(max(adf))
    			fmt[0] = f"{gt}/{gt}"
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
