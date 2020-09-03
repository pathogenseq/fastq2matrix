import sys
import fastq2matrix as fm
import argparse
import csv
import json
from collections import defaultdict

drugs = ['rifampicin', 'isoniazid', 'ethambutol', 'pyrazinamide', 'streptomycin', 'amikacin', 'kanamycin', 'capreomycin', 'ofloxacin', 'moxifloxacin', 'ciprofloxacin', 'levofloxacin', 'cycloserine', 'ethionamide', 'para-aminosalicylic_acid', 'rifabutin', 'bedaquiline', 'delamanid', 'clofazimine','linezolid','clarithromycin']

def main(args):
    samples = []
    for l in open(args.samples):
        samples = [x.rstrip() for x in open(args.samples).readlines()]

    for s in samples:
        fm.filecheck("per_sample/%s%s" % (s,args.alignment_extension))

    if fm.nofolder("%(dir)s/kraken" % vars(args)):
        fm.run_cmd("%(dir)s/kraken" % vars(args))

    args.cmd_file = fm.get_random_file()
    with open(args.cmd_file,"w") as O:
        for s in samples:
            args.sample = s

            if fm.nofile("%(dir)s/per_sample/%(sample)s.median_dp.txt" % vars(args)):
                O.write("printf %s\"\\t\"$(bedtools genomecov -d -ibam %s/per_sample/%s%s | datamash median 3)\"\\n\" > %s/per_sample/%s.median_dp.txt\n" % (s,args.dir,s,args.alignment_extension,args.dir,s))

            if fm.nofile("%(dir)s/kraken/%(sample)s.done" % vars(args)):
                O.write("kraken2 --db /run/user/506/standard --gzip-compressed --paired %(dir)s/fastq/%(sample)s_1.fastq.gz %(dir)s/fastq/%(sample)s_2.fastq.gz --report %(dir)s/kraken/%(sample)s.report.txt --out %(dir)s/kraken/%(sample)s.out.txt --threads 10 --memory-mapping && touch %(dir)s/kraken/%(sample)s.done\n" % vars(args))

    fm.run_cmd("cat %(cmd_file)s | parallel -j %(io_heavy_threads)s" % vars(args), verbose=2)
    fm.rm_files([args.cmd_file])

    sample_metrics = []
    for s in samples:
        res = {"sample":s}
        args.sample = s
        for i,l in enumerate(open("%(dir)s/per_sample/%(sample)s.bqsr.bamstats" % vars(args))):
            row = l.rstrip().split()
            if i in [2,3]:
                res[row[3]] = int(row[0])
            elif i==4:
                res[row[3]] = int(row[0])
                res["mapped_percent"] = float(row[4].replace("(","").replace("%",""))
            else:
                pass

        kraken_results = {}
        for l in open("%(dir)s/kraken/%(sample)s.report.txt" % vars(args)):
            row = l.strip().split()
            if row[3] not in kraken_results:
                kraken_results[row[3]] = (float(row[0])," ".join(row[5:]))
            if float(row[0]) > kraken_results[row[3]][0]:
                kraken_results[row[3]] = (float(row[0])," ".join(row[5:]))
        res["kraken_genus"] = "%s (%.2f)" % (kraken_results["G"][1],kraken_results["G"][0])
        res["kraken_genus1"] = "%s (%.2f)" % (kraken_results["G1"][1],kraken_results["G1"][0])
        res["kraken_species"] = "%s (%.2f)" % (kraken_results["S"][1],kraken_results["S"][0])

        tbprofiler_result = json.load(open("%(dir)s/tbprofiler/results/%(sample)s.results.json" % vars(args)))
        res["lineage"] = tbprofiler_result["main_lin"]
        res["sub-lineage"] = tbprofiler_result["sublin"]
        res["drtype"] = tbprofiler_result["drtype"]
        tmp_drugs = defaultdict(list)
        for var in tbprofiler_result["dr_variants"]:
            for d in var["drugs"]:
                tmp_drugs[d["drug"]].append("%s_%s (%.2f)" % (var["gene"],var["change"],var["freq"]))
        for d in drugs:
            res[d] = ", ".join(tmp_drugs[d])

        sample_metrics.append(res)

    with open(args.out+".sample_info.csv","w") as O:
        writer = csv.DictWriter(O, fieldnames=list(sample_metrics[0]))
        writer.writeheader()
        writer.writerows(sample_metrics)

    vcf = fm.vcf_class(args.vcf)
    if fm.nofile(args.vcf+".stats.txt"):
        fm.run_cmd("bcftools norm -m - -f %(ref)s %(vcf)s | bcftools stats -v -s - > %(vcf)s.stats.txt" % (vars(args)))

    vcf_stats = vcf.load_stats()

    results = {
        "number of samples": vcf_stats["number of samples"],
        "number of records": vcf_stats["number of records"],
        "number of SNPs": vcf_stats["number of SNPs"],
        "number of indels": vcf_stats["number of indels"],
    }

    snp_results = []
    if fm.nofile(args.vcf+".csq_info.txt"):
        fm.run_cmd("bcftools view -V indels %(vcf)s | bcftools norm -m - -f %(ref)s | bcftools csq -f %(ref)s -g %(gff)s | correct_tb_csq.py | bcftools +fill-tags | bcftools query -f '%%POS\\t%%REF\\t%%ALT\\t%%AF\\t%%AC\\t%%BCSQ\\n' > %(vcf)s.csq_info.txt" % vars(args))
        fm.run_cmd("bcftools view -v indels %(vcf)s | bcftools norm -m - -f %(ref)s | bcftools csq -f %(ref)s -g %(gff)s | correct_tb_csq.py | bcftools +fill-tags | bcftools query -f '%%POS\\t%%REF\\t%%ALT\\t%%AF\\t%%AC\\t%%BCSQ\\n' >> %(vcf)s.csq_info.txt" % vars(args))

    variant_info = vcf.get_variant_data(args.ref,args.gff)
    with open(args.out+".variant_info.csv","w") as O:
        writer = csv.DictWriter(O, fieldnames=list(variant_info[0]))
        writer.writeheader()
        writer.writerows(variant_info)


parser = argparse.ArgumentParser(description='XXX pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--out',help='Root directory file',required=True)
parser.add_argument('--vcf',help='VCF',required=True)
parser.add_argument('--ref',help='Reference sequence',required=True)
parser.add_argument('--gff',help='Reference sequence',required=True)
parser.add_argument('--dir',default=".",help='Root directory file')
parser.add_argument('--samples',help='Root directory file',required=True)
parser.add_argument('--alignment-extension',default=".bqsr.bam",help='Root directory file')
parser.add_argument('--io-heavy-threads',default=4,help='Root directory file')
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
