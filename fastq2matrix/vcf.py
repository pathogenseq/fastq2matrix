import os
from .utils import run_cmd, nofile, get_random_file, add_arguments_to_self, index_vcf
from collections import defaultdict
import uuid

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
        index_vcf(filename,threads=threads)
        self.temp_file = get_random_file()
        run_cmd("bcftools query -l %(filename)s > %(temp_file)s" % vars(self))
        for l in open(self.temp_file):
            self.samples.append(l.rstrip())
        os.remove(self.temp_file)

    def vcf_to_fasta(self,ref_file,threads=4,chunk_size = 50000, nofilt=True):
        self.ref_file = ref_file
        self.window = chunk_size-1
        self.step = chunk_size
        self.cmd_split_chr = "bedtools makewindows -g %(ref_file)s.fai -w %(window)s -s %(step)s | awk '{print $1\":\"$2\"-\"$3}'" % vars(self)
        self.tmp_file = "%s.tmp.txt" % self.prefix
        self.threads = threads
        if nofilt:
            cmd = "%(cmd_split_chr)s | parallel -j %(threads)s \"bcftools view  %(filename)s -r {} | bcftools view -V indels | bcftools query -f '%%POS[\\t%%IUPACGT]\\n' | sed 's/\.[\/|]\./N/g' |  datamash transpose > %(prefix)s.{}.tmp.txt\"" % vars(self)
        else:
            cmd = "%(cmd_split_chr)s | parallel -j %(threads)s \"bcftools view  %(filename)s -r {} | bcftools view -V indels | setGT.py | bcftools view -a | bcftools filter -e 'GT=\\\"het\\\"' -S . | bcftools view -i 'F_PASS(GT!=\\\"mis\\\")>0.9' | bcftools view -c 1 | bcftools +fill-tags | bcftools view -e 'AF==1 || AF==0' | bcftools norm -f %(ref_file)s | bcftools query -f '%%POS[\\t%%IUPACGT]\\n' | sed 's/\.[\/|]\./N/g' |  datamash transpose > %(prefix)s.{}.tmp.txt\"" % vars(self)
        run_cmd(cmd)
        cmd = "paste `%(cmd_split_chr)s | awk '{print \"%(prefix)s.\"$1\".tmp.txt\"}'` > %(tmp_file)s" % vars(self)
        run_cmd(cmd)
        cmd = "rm `%(cmd_split_chr)s | awk '{print \"%(prefix)s.\"$1\".tmp.txt\"}'`" % vars(self)
        run_cmd(cmd)
        with open(self.prefix+".snps.fa","w") as O:
            for i,l in enumerate(open(self.tmp_file)):
                row = l.rstrip().split()
                if i==0: continue
                s = self.samples[i-1]
                seq = "".join(row)
                O.write(">%s\n%s\n" % ( s,seq))
        run_cmd("rm %s" % self.tmp_file)
        return self.prefix+".snps.fa"

    def vcf_to_matrix(self, iupacgt=True, na = "N"):
        self.matrix_file = self.prefix+".mat"
        self.binary_matrix_file = self.prefix+".mat.bin"
        self.na = na
        if iupacgt:
            self.matrix_file = self.prefix+".noniupac.mat"
            O = open(self.matrix_file,"w").write("chr\tpos\tref\t%s\n" % ("\t".join(self.samples)))
            run_cmd("bcftools query -f '%%CHROM\\t%%POS\\t%%REF[\\t%%TGT]\\n' %(filename)s | sed 's/\.\/./%(na)s/g; s/\([ACTG]\)\///g; s/|//g' | sed -r 's/([ACGT])\\1+/\\1/g' >> %(matrix_file)s" % vars(self))
        else:
            O = open(self.matrix_file,"w").write("chr\tpos\tref\t%s\n" % ("\t".join(self.samples)))
            run_cmd("bcftools query -f '%%CHROM\\t%%POS\\t%%REF[\\t%%IUPACGT]\\n' %(filename)s | tr '|' '/' | sed 's/\.\/\./%(na)s/g' >> %(matrix_file)s" % vars(self))

        O = open(self.binary_matrix_file,"w").write("chr\tpos\tref\t%s\n" % ("\t".join(self.samples)))
        run_cmd("bcftools query -f '%%CHROM\\t%%POS\\t%%REF[\\t%%GT]\\n' %(filename)s | tr '|' '/' | sed 's/\.\/\./%(na)s/g' | sed 's/0\/1/0.5/g' | sed 's/[123456789]\/[123456789]/1/g' | sed 's/0\/0/0/g' >> %(binary_matrix_file)s" % vars(self))

    def get_plink_dist(self,pca=True,mds=True):
        self.tempfile = get_random_file(extension=".vcf")
        run_cmd("bcftools view %(filename)s > %(tempfile)s" % vars(self))
        run_cmd("plink --vcf %(tempfile)s --distance square --double-id --allow-extra-chr --out %(prefix)s.temp" % vars(self))
        O = open("%(prefix)s.dist" % vars(self),"w")
        dists = []
        for l in open("%(prefix)s.temp.dist" % vars(self)):
            row = [float(d)/2 for d in l.rstrip().split()]
            O.write("%s\n" % "\t".join([str(x) for x in row]))
            dists.append(row)
        O.close()
        if pca:
            run_cmd("plink --vcf %(tempfile)s --pca --double-id --allow-extra-chr --out %(prefix)s.pca" % vars(self))
        if mds:
            run_cmd("plink --vcf %(tempfile)s --mds-plot 10 eigendecomp --cluster --double-id --allow-extra-chr --out %(prefix)s.pca" % vars(self))

        run_cmd("rm %(tempfile)s* %(prefix)s.temp* %(prefix)s.pca.log %(prefix)s.pca.nosex" % vars(self))
        return dists

    def load_stats(self,convert=False,ref=None):
            add_arguments_to_self(self,locals())
            self.smallest_bin = 1/len(self.samples)
            self.stats_file = "%s.stats.txt" % self.filename
            if nofile(self.stats_file):
                run_cmd("bcftools stats -v -s - %(filename)s > %(stats_file)s" % vars(self))

            results = defaultdict(lambda:defaultdict(dict))
            for l in open(self.stats_file):
                row = l.rstrip().split("\t")
                if l[0]=="#": continue
                if row[0]=="SN":
                    results["SN"][row[2][:-1]] = int(row[3])
                elif row[0]=="AF":
                    results["AF"]["SNP"][float(row[2])] = int(row[3])
                    results["AF"]["INDEL"][float(row[2])] = int(row[6])
                elif row[0]=="QUAL":
                    results["QUAL"]["SNP"][int(row[2])] = int(row[3])
                    results["QUAL"]["INDEL"][int(row[2])] = int(row[6])
                elif row[0]=="IDD":
                    results["IDD"][int(row[2])] = int(row[3])
                elif row[0]=="ST":
                    results["ST"][row[2]] = int(row[3])
                elif row[0]=="DP":
                    if row[2][0]==">": continue
                    results["DP"][int(row[2])] = int(row[3])
                elif row[0]=="PSC":
                    results["PSC"][row[2]]["nRefHom"] = int(row[3])
                    results["PSC"][row[2]]["nNonRefHom"] = int(row[4])
                    results["PSC"][row[2]]["nHets"] = int(row[5])
            return results

    def get_variant_data(self,ref,gff):
        add_arguments_to_self(self,locals())

        if nofile(self.filename+".csq_info.txt"):
            run_cmd("bcftools view -V indels %(filename)s | bcftools norm -m - -f %(ref)s | bcftools csq -f %(ref)s -g %(gff)s | correct_tb_csq.py | bcftools +fill-tags | bcftools query -f '%%POS\\t%%REF\\t%%ALT\\t%%AF\\t%%AC\\t%%BCSQ\\n' > %(filename)s.csq_info.txt" % vars(self))
            run_cmd("bcftools view -v indels %(filename)s | bcftools norm -m - -f %(ref)s | bcftools csq -f %(ref)s -g %(gff)s | correct_tb_csq.py | bcftools +fill-tags | bcftools query -f '%%POS\\t%%REF\\t%%ALT\\t%%AF\\t%%AC\\t%%BCSQ\\n' >> %(filename)s.csq_info.txt" % vars(self))

        results = defaultdict(list)
        for l in open("%(filename)s.csq_info.txt" % vars(self)):
            print(l)
            pos,ref,alt,af,ac,csq_str = l.strip().split()
            alt_af = float(af)
            alt_ac = int(int(ac)/2)
            csqs = csq_str.split(",")
            types = set()
            changes = set()
            genes = set()
            nucleotide_changes = set()
            pos = int(pos)
            tmp = {"pos":pos, "ref":ref, "alts":alt, "alt_af":alt_af, "alt_ac":alt_ac}
            print(l)
            for i in range(len(csqs)):
                if csqs[i][0]=="@":
                    tmp["csq_link"] = csqs[i]
                elif csqs[i]==".":
                    types.add("non_coding")
                    changes.add("NA")
                    genes.add("NA")
                    nucleotide_changes.add("NA")
                else:
                    csq = csqs[i].split("|")
                    types.add(csq[0].replace("*",""))
                    changes.add(csq[5])
                    genes.add(csq[1])
                    nucleotide_changes.add(csq[6])

            tmp["types"] = ",".join(types)
            tmp["changes"] = ",".join(changes)
            tmp["gene"] = ",".join(genes)
            tmp["nucleotide_change"] = ",".join(nucleotide_changes)
            results[pos].append(tmp)

        variant_rows = []
        for pos in results:
            for row in results[pos]:
                if "csq_link" in row:
                    tmp = [x for x in results[int(row["csq_link"][1:])] if str(pos) in str(x)]
                    if len(tmp)!=1:
                        import pdb; pdb.set_trace()
                    else:
                        if len(results[int(row["csq_link"][1:])])==1:
                            variant_rows.append(tmp[0].copy())
                        t = tmp[0]
                        t["pos"] = row["pos"]
                        t["alts"] = row["alts"]
                        t["alt_af"] = row["alt_af"]
                        t["alt_ac"] = row["alt_ac"]
                        # variant_rows.append(t)
                else:
                    variant_rows.append(row)
        # import pdb; pdb.set_trace()

        return sorted(variant_rows,key=lambda x:x["pos"])
    def filter_by_af(self,maf,pop_file=None,threads=4):
        add_arguments_to_self(self,locals())

        tmp_pop_file = self.uuid+".pops"
        with open(tmp_pop_file,"w") as O:
            for l in open(pop_file):
                row = l.strip().split()
                O.write("%s\t%s\t%s\n" % (row[0],row[0],row[1]))


        run_cmd("bcftools norm --threads %(threads)s -m - %(filename)s  | bcftools annotate --threads %(threads)s --set-id '%%POS' -Oz -o %(uuid)s.vcf.gz" % vars(self))
        # run_cmd("plink --vcf %(uuid)s.vcf.gz --make-bed --allow-extra-chr --keep-allele-order --out %(uuid)s" % vars(self))
        self.within = "--within %s" % tmp_pop_file if pop_file else ""
        run_cmd("plink --vcf %(uuid)s.vcf.gz --freq %(within)s  --allow-extra-chr --out %(uuid)s" % vars(self))

        with open(self.uuid+".extract.bed","w") as O:
            for l in open(self.uuid+".frq.strat"):
                row = l.strip().split()
                if row[0]=="CHR": continue
                if float(row[5])>self.maf:
                        O.write("%s\t%s\t%s\n" % (row[0],int(row[1])-1,row[1]))

        run_cmd("bcftools view --threads %(threads)s -T %(uuid)s.extract.bed %(filename)s -Oz -o %(prefix)s.pop_maf_filt_%(maf)s.vcf.gz" % vars(self))
        run_cmd("mv %(uuid)s.frq.strat %(prefix)s.frq.strat" % vars(self))
        run_cmd("rm %(uuid)s*" % vars(self))
