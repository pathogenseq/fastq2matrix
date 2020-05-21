from fastq2matrix import run_cmd, nofile, nofolder, vcf_class, get_random_file
import argparse



def main(args):
    vcf_obj = vcf_class(args.vcf)
    if args.snps:
        vcf_obj.vcf_to_fasta(args.ref,nofilt=args.snps_no_filt)
    else:
        args.sample_file = get_random_file()
        open(args.sample_file,"w").write("\n".join(vcf_obj.samples)+"\n")
        run_cmd('cat %(sample_file)s | parallel  --bar -j %(threads)s "bcftools consensus -f %(ref)s -s {} %(vcf)s | sed \'s/^>.*/>{}/\' > {}.tmp.fasta"' % vars(args))
        run_cmd('cat %s > %s.fa' % (" ".join(["%s.tmp.fasta" % s for s in vcf_obj.samples]), vcf_obj.prefix))
        run_cmd('rm %s %s' % (" ".join(["%s.tmp.fasta" % s for s in vcf_obj.samples]), args.sample_file))



parser = argparse.ArgumentParser(description='VCF mergin pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--vcf',help='sample file',required=True)
parser.add_argument('--ref',help='reference file',required=True)
parser.add_argument('--snps',action="store_true",help='Only use SNPs')
parser.add_argument('--snps-no-filt',action="store_true",help='Only use SNPs')
parser.add_argument('--threads','-t',default=4,help='Number of threads')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
