#!/usr/bin/python3

import sys
import os

errors = 0
def check_log(fn):
    global errors
    for line in open(fn):
        if line.startswith("Finished"):
            print("ok %s" % fn)
            return
    errors += 1
    print("ERROR %s Finished not found" % fn)

d = "../test_logs/"
fns = os.listdir(d)
for fn in fns:
    if fn.endswith(".log"):
        check_log(d + fn)

exit(1 if errors > 0 else 0)
