from fastq2matrix import run_cmd, nofile, nofolder, vcf_class, get_random_file
import argparse



def main(args):
    vcf_obj = vcf_class(args.vcf)
    if args.whole_genome:
        args.sample_file = get_random_file()
        open(args.sample_file,"w").write("\n".join(vcf_obj.samples)+"\n")
        run_cmd('cat %(sample_file)s | parallel  --bar -j %(threads)s "bcftools consensus -f %(ref)s -s {} %(vcf)s | sed \'s/^>.*/>{}/\' > {}.tmp.fasta"' % vars(args))
        run_cmd('cat %s > %s.fa' % (" ".join(["%s.tmp.fasta" % s for s in vcf_obj.samples]), vcf_obj.prefix))
        run_cmd('rm %s %s' % (" ".join(["%s.tmp.fasta" % s for s in vcf_obj.samples]), args.sample_file))
        if args.tree:
            run_cmd("iqtree -s %s.fa -m GTR+G+ASC -nt AUTO" % vcf_obj.prefix)
    else:
        fasta_file = vcf_obj.vcf_to_fasta(args.ref,nofilt=args.snps_no_filt)
        if args.tree:
            run_cmd("iqtree -s %s -m GTR+G+ASC -nt AUTO" % fasta_file)




parser = argparse.ArgumentParser(description='VCF mergin pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--vcf',help='sample file',required=True)
parser.add_argument('--ref',help='reference file',required=True)
parser.add_argument('--snps',action="store_true",help='Only use SNPs (redundant)')
parser.add_argument('--tree',action="store_true",help='Generate phylogenetic tree')
parser.add_argument('--whole-genome',action="store_true",help='Generate whole genome sequences')
parser.add_argument('--snps-no-filt',action="store_true",help='Do not filter snps')
parser.add_argument('--threads','-t',default=4,help='Number of threads')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
