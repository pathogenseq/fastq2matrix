import sys
import argparse
import ete3
import subprocess
from collections import defaultdict
from tqdm import tqdm
import fastq2matrix as fm


def main(args):

    vcf_class = fm.vcf(args.vcf)
    vcf_positions = vcf_class.get_positions()

    if not args.fasta:
        if not args.ref:
            sys.stderr.write("\nERROR: Please supply a reference with --ref\n\n")
            quit()
        fm.run_cmd("vcf2fasta.py --vcf %(vcf)s --snps --ref %(ref)s --snps-no-filt" % vars(args))
        args.fasta = "%s.snps.fa" % vcf_class.prefix
    if fm.nofile("%s.asr.state" % args.fasta):
        fm.run_cmd("iqtree -m %(model)s -te %(tree)s -s %(fasta)s -nt AUTO -asr -pre %(fasta)s.asr" % vars(args))

    tree = ete3.Tree("%s.asr.treefile" % args.fasta,format=1)
    node_names = set([tree.name] + [n.name.split("/")[0] for n in tree.get_descendants()])
    leaf_names = set(tree.get_leaf_names())
    internal_node_names = node_names - leaf_names

    states_file = "%s.asr.state" % args.fasta
    states = defaultdict(dict)
    sys.stderr.write("Loading states\n")
    for l in tqdm(open(states_file)):
        if l[0] == "#": continue
        row = l.strip().split()
        if row[0] == "Node": continue
        site = int(row[1])
        if row[0] not in internal_node_names: continue
        states[site][row[0]] = row[2]

    seqs = fm.fasta(args.fasta).fa_dict
    for site in tqdm(list(states)):
        for sample in seqs:
            states[site][sample] = seqs[sample][site-1]


    acgt = set(["A", "C", "G", "T", "a", "c", "g", "t"])
    convergent_sites = []
    for site in tqdm(list(states)):
        nucleotides = set([states[site][n] for n in node_names])
        if len(nucleotides)==1: continue

        # Set up storage objects
        origins = []

        tree.add_feature("state", states[site][tree.name])
        for n in tree.traverse():
            if n == tree: continue
            node_state = states[site][n.name]
            if node_state != n.get_ancestors()[0].state and node_state in acgt and n.get_ancestors()[0].state in acgt:
                origins.append(n.name)
            n.add_feature("state", node_state)
        if len(origins) > 1:
            convergent_sites.append((site, vcf_positions[site-1], origins))

    with open(args.out,"w") as O:
        for site in convergent_sites:
            O.write("%s\t%s\n" % (site[1][1],len(site[2])))




parser = argparse.ArgumentParser(description='Assembly to VCF',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--vcf',help='VCF file with just SNPS',required=True)
parser.add_argument('--tree',help='Tree file',required=True)
parser.add_argument('--out',help='Output file',required=True)
group = parser.add_mutually_exclusive_group()
group.add_argument('--fasta',help='Fasta file contining exactly the same number of SNPs as in the VCF. If you want to be safe that you have exactly the same just use the --ref command instead to generate a SNPs fasta')
group.add_argument('--ref',help='Reference file. Use this if you want to create a SNPs fasta from the VCF')
parser.add_argument('--model',default="GTR+G",help='Model used by iqtree')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
