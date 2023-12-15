#!/usr/bin/python
import sys, os, glob, re
group = sys.argv[1]
path = "/ccs/home/arodriguez/med112/task0101113/output/HARE_ANC_Run"
log_files = glob.glob("%s/*/%s/step1/*.%s.out.gpu.log" % (path, group, group))
step1_times = {}
for i in log_files:
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
        step1_times[i] = phecode_times[-1]

for i in step1_times:
    print("%s\t%.2f\tminutes" % (i.split("/")[8], float(step1_times[i])/60))
