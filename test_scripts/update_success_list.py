#!/usr/bin/python3

import sys

ver = sys.argv[1].strip()
date = sys.argv[2].strip()

fn = "../test_results/success_list.txt"
lines = []
exists = False
try:
    for line in open(fn):
        lines.append(line[:-1])
    exists = True
except:
    pass

found = False
for line in lines:
    flds = line.split('\t')
    if flds[0] == ver:
        print("Previously suceeded %s" % line)
        sys.exit(0)

f = open(fn, "a" if exists else "w")
line = ver + "\t" + date
f.write(line + "\n")
f.close()
print("Success %s" % line)
