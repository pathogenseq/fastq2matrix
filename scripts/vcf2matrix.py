#! /usr/bin/env python
import sys
import subprocess
import argparse
import random
import os
rand_generator = random.SystemRandom()

def get_random_file(prefix = None,extension=None):
	randint = rand_generator.randint(1,999999)
	if prefix:
		if extension:
			return "%s.%s.%s" % (prefix,randint,extension)
		else:
			return "%s.%s.txt" % (prefix,randint)
	else:
		if extension:
			return "%s.tmp.%s" % (randint,extension)
		else:
			return "%s.tmp.txt" % (randint)

def log(msg,ext=False):
	sys.stderr.write("\n"+str(msg)+"\n")
	if ext:
		exit(1)

def run_cmd(cmd,verbose=1,target=None):
	"""
	Wrapper to run a command using subprocess with 3 levels of verbosity and automatic exiting if command failed
	"""
	if target and filecheck(target): return True
	cmd = "set -u pipefail; " + cmd
	if verbose==2:
		sys.stderr.write("\nRunning command:\n%s\n" % cmd)
		stdout = open("/dev/stdout","w")
		stderr = open("/dev/stderr","w")
	elif verbose==1:
		sys.stderr.write("\nRunning command:\n%s\n" % cmd)
		stdout = open("/dev/null","w")
		stderr = open("/dev/null","w")
	else:
		stdout = open("/dev/null","w")
		stderr = open("/dev/null","w")

	res = subprocess.call(cmd,shell=True,stderr = stderr,stdout = stdout)
	stderr.close()
	if res!=0:
		print("Command Failed! Please Check!")
		exit(1)

def nofile(filename):
	"""
	Return True if file does not exist
	"""
	if not os.path.isfile(filename):
		return True
	else:
		return False

def get_vcf_prefix(filename):
	if filename[-4:]==".bcf":
		return filename[:-4]
	elif filename[-5:]==".gbcf":
		return filename[:-5]
	elif filename[-7:]==".vcf.gz":
		return filename[:-7]
	elif filename[-4:]==".vcf":
		return filename[:-4]
	else:
		return filename

class vcf_class:
	def __init__(self,filename,threads=4):
		self.samples = []
		self.filename = filename
		self.threads = threads
		self.prefix = get_vcf_prefix(filename)
		if nofile(filename+".csi"):
			run_cmd("bcftools index  %(filename)s" % vars(self))
		self.temp_file = get_random_file()
		run_cmd("bcftools query -l %(filename)s > %(temp_file)s" % vars(self))
		for l in open(self.temp_file):
			self.samples.append(l.rstrip())
		os.remove(self.temp_file)
	def vcf_to_matrix(self,):
		self.matrix_file = self.prefix+".mat"
		self.binary_matrix_file = self.prefix+".mat.bin"
		O = open(self.matrix_file,"w").write("chr\tpos\tref\t%s\n" % ("\t".join(self.samples)))
		run_cmd("bcftools query -f '%%CHROM\\t%%POS\\t%%REF[\\t%%IUPACGT]\\n' %(filename)s | tr '|' '/' | sed 's/\.\/\./N/g' >> %(matrix_file)s" % vars(self))
		O = open(self.binary_matrix_file,"w").write("chr\tpos\tref\t%s\n" % ("\t".join(self.samples)))
		run_cmd("bcftools query -f '%%CHROM\\t%%POS\\t%%REF[\\t%%GT]\\n' %(filename)s | tr '|' '/' | sed 's/\.\/\./N/g' | sed 's/0\/1/0.5/g' | sed 's/1\/1/1/g' | sed 's/0\/0/0/g' >> %(binary_matrix_file)s" % vars(self))

def main(args):
	if nofile(args.vcf): quit("Can't find %s... Exiting!" % args.vcf)
	vcf = vcf_class(args.vcf)
	vcf.vcf_to_matrix()

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--vcf',help='VCF file',required=True)
parser.add_argument('--threads',default=4, type=int, help='Number of threads for parallel operations')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
