import sys
import re

for l in sys.stdin:
    row = l.rstrip().split("\t")
    if "BCSQ=" in l:
        info = row[7].split(";")
        for i in range(len(info)):
            if info[i][:4]=="BCSQ":
                new_bcsq = []
                for bcsq in info[i].split(","):
                    if "@" in bcsq:
                        new_bcsq.append(bcsq)
                    else:
                        tmp = bcsq.split("|")
                        if len(tmp)==7:
                            new_bcsq.append("|".join(tmp))
                        else:
                            for x in range(7-len(tmp)):
                                tmp = tmp + ["NA","NA","NA"]
                            new_bcsq.append("|".join(tmp))
                info[i] = ",".join(new_bcsq)
        row[7] = ";".join(info)
    sys.stdout.write("\t".join(row)+"\n")
