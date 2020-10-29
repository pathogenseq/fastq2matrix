import sys
import gzip
import os.path
import subprocess
import csv
from collections import defaultdict
import json
import random
import math
import re
import uuid
rand_generator = random.SystemRandom()

def chunk_reference(ref,n):
    genome_chunks = []
    for l in cmd_out(f"bedtools makewindows -n {n} -g {ref}.fai"):
        row = l.strip().split()
        genome_chunks.append("%s:%s-%s" % (row[0],row[1],row[2]))
    return genome_chunks

def get_contigs_from_fai(filename):
    contigs = []
    for l in open(filename):
        contigs.append(l.strip().split()[0])
    return contigs

def debug(x):
    sys.stderr.write("#"*40+"\n")
    sys.stderr.write(x+"\n")
    sys.stderr.write("#"*40+"\n")

def filetype(x):
    for l in cmd_out("file %s" % x):
        pass
    row = l.rstrip().split()
    return row[1]

def revcom(s):
        """Return reverse complement of a sequence"""
        def complement(s):
                        basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
                        letters = list(s)
                        letters = [basecomplement[base] for base in letters]
                        return ''.join(letters)
        return complement(s[::-1])

def stdev(arr):
    mean = sum(arr)/len(arr)
    return math.sqrt(sum([(x-mean)**2 for x in arr])/len(arr))

def add_arguments_to_self(self,args):
    for x in args:
        if x=="self": continue
        vars(self)[x] = args[x]
    if "kwargs" in args:
        for x in args["kwargs"]:
            vars(self)[x] = args["kwargs"][x]
    vars(self)["uuid"] = str(uuid.uuid4())

def cmd_out(cmd,verbose=1):
    cmd = "set -u pipefail; " + cmd
    if verbose==2:
        sys.stderr.write("\nRunning command:\n%s\n" % cmd)
        stderr = open("/dev/stderr","w")
    elif verbose==1:
        sys.stderr.write("\nRunning command:\n%s\n" % cmd)
        stderr = open("/dev/null","w")
    else:
        stderr = open("/dev/null","w")
    try:
        res = subprocess.Popen(cmd,shell=True,stderr = stderr,stdout=subprocess.PIPE)
        for l in res.stdout:
            yield l.decode().rstrip()
    except:
        sys.stderr.write("Command Failed! Please Check!")
        exit(1)
    stderr.close()

def get_random_file(prefix = None,extension=None):
    randint = rand_generator.randint(1,999999)
    if prefix:
        if extension:
            return "%s.%s%s" % (prefix,randint,extension)
        else:
            return "%s.%s.txt" % (prefix,randint)
    else:
        if extension:
            return "%s.tmp%s" % (randint,extension)
        else:
            return "%s.tmp.txt" % (randint)

def log(msg,ext=False):
    sys.stderr.write("\n"+str(msg)+"\n")
    if ext:
        exit(1)

def init_params():
    conf = json.load(open("%s/%s" % (sys.prefix,"pathogenseq.conf")))
    return conf

def filecheck(filename):
    """
    Check if file is there and quit if it isn't
    """
    if not os.path.isfile(filename):
        sys.stderr.write("Can't find %s\n" % filename)
        exit(1)
    else:
        return filename

def foldercheck(filename):
    """
    Check if file is there and quit if it isn't
    """
    if not os.path.isdir(filename):
        sys.stderr.write("Can't find %s\n" % filename)
        exit(1)
    else:
        return filename

def debug(s):
    sys.stderr.write("#"*40+"\n")
    sys.stderr.write("%s\n" % s)
    sys.stderr.write("#"*40+"\n")

def nofile(filename):
    """
    Return True if file does not exist
    """
    if not os.path.isfile(filename):
        return True
    else:
        return False

def nofolder(filename):
    """
    Return True if file does not exist
    """
    if not os.path.isdir(filename):
        return True
    else:
        return False

def create_seq_dict(ref):
    if nofile("%s.dict" % (ref.replace(".fasta","").replace(".fa",""))):
        run_cmd("gatk CreateSequenceDictionary -R %s" % ref)

def faidx(ref):
    if nofile("%s.fai" % (ref)):
        run_cmd("samtools faidx %s" % ref)

def tabix_vcf(vcf):
    if nofile("%s.tbi" % (vcf)):
        run_cmd("bcftools index -t %s" % vcf)

def bwa_index(ref):
    """
    Create BWA index for a reference
    """
    if nofile("%s.bwt"%ref):
        cmd = "bwa index %s" % ref
        run_cmd(cmd)

def run_cmd(cmd,verbose=1,target=None):
    """
    Wrapper to run a command using subprocess with 3 levels of verbosity and automatic exiting if command failed
    """
    if target and     (target): return True
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
        sys.stderr.write("Command Failed! Please Check!")
        exit(1)

def index_bam(bamfile,threads=4,overwrite=False):
    """
    Indexing a bam file
    """
    cmd = "samtools index -@ %s %s" % (threads,bamfile)
    if     (bamfile):
        if nofile(bamfile+".bai"):
            run_cmd(cmd)
        elif os.path.getmtime(bamfile+".bai")<os.path.getmtime(bamfile) or overwrite:
            run_cmd(cmd)



def index_vcf(bcffile,threads=4,overwrite=False):
    """
    Indexing a bam file
    """
    cmd = "bcftools index --threads %s -f %s" % (threads,bcffile)
    if     (bcffile):
        if nofile(bcffile+".csi"):
            run_cmd(cmd)
        elif os.path.getmtime(bcffile+".csi")<os.path.getmtime(bcffile) or overwrite:
            run_cmd(cmd)

def verify_fq(filename):
    """
    Return True if input is a valid fastQ file
    """
    FQ = open(filename) if filename[-3:]!=".gz" else gzip.open(filename)
    l1 = FQ.readline()
    if l1[0]!="@":
        sys.stderr.write("First character is not \"@\"\nPlease make sure this is fastq format\nExiting...")
        exit(1)
    else:
        return True

def rm_files(x,verbose=True):
    """
    Remove a files in a list format
    """
    for f in x:
        if os.path.isfile(f):
            if verbose: sys.stderr.write("Removing %s\n" % f)
            os.remove(f)

def download_from_ena(acc):
    if len(acc)==9:
        dir1 = acc[:6]
        cmd = "wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/%s/%s*" % (dir1,acc,acc)
    elif len(acc)==10:
        dir1 = acc[:6]
        dir2 = "00"+acc[-1]
        cmd = "wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/%s/%s/%s/%s*" % (dir1,dir2,acc,acc)
    else:
        sys.stderr.write("Check Accession: %s" % acc)
        exit(1)
    run_cmd(cmd)
