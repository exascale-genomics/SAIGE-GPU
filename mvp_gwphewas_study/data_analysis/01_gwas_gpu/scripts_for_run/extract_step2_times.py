#!/usr/bin/python
import sys, os, glob, re
group = sys.argv[1]
path = "./output/HARE_ANC_Run"
phe_dirs = glob.glob("%s/*/%s/step2" % (path, group))
for j in phe_dirs:
    step2_times = {}
    log_files = glob.glob("%s/*.%s.assoc.gpu.txt.log" % (j, group))
    qty = len(log_files)
    for i in log_files:
        #print(i)
        fh = open(i, "r")
        phecode_times = []
        for line in fh:
            if "elapsed" in line:
                line = fh.readline().strip()
                #print(line)
                if not re.search('[a-zA-Z\[\]:]', line) and len(re.split(' +', line.rstrip("\n"))) == 3:
                    phecode_times.append(re.split(' +', line.rstrip("\n"))[2])
        #print(phecode_times)
        if len(phecode_times) > 0:
            phecode_times.sort(key=lambda x: (float(x), len(x)))
            step2_times[i] = phecode_times[-1]
        fh.close()

    for i in step2_times:
        print("%s\t%.2f\tminutes\t%s" % (i.split("/")[8], float(step2_times[i])/60, qty))
