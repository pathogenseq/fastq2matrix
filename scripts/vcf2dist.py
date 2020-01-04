from fastq2matrix import vcf_class, run_cmd, nofile, nofolder, vcf_class, get_random_file
import argparse

def main(args):
	vcf_obj = vcf_class(args.vcf)
	run_cmd("plink --vcf %(vcf)s --distance square --double-id --allow-extra-chr --vcf-half-call missing --out %(vcf)s" % vars(args))
	outfile = open("%s.dists" % vcf_obj.prefix,"w")
	outfile.write("%s\n" % "\t".join(vcf_obj.samples))
	for l in open("%s.dist" % args.vcf):
		row = l.strip().split()
		outfile.write("%s\n" % "\t".join([str(float(x)/2) for x in row]))
	outfile.close()
	run_cmd("rm %(vcf)s.dist %(vcf)s.log" % vars(args))




parser = argparse.ArgumentParser(description='VCF mergin pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--vcf',help='sample file',required=True)
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
