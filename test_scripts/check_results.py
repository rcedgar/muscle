#!/usr/bin/python3

import sys
import re

errors = 0
def qscore(fn, Q=None, TC=None):
	global errors
	q = None
	tc = None
	lines = []
	try:
		for line in open(fn):
			lines.append(line)

	except:
		print("ERROR reading file " + fn)
		errors += 1
		return
	for line in lines:
		if line.find(", TC=") > 0:
			M = re.search(r"Q=([0-9.]*), TC=([0-9.]*)", line)
			if not M is None:
				q = float(M.group(1))
				tc = float(M.group(2))
				break
	if q is None or tc is None:
		print("ERROR Q= TC= not found " + fn)
		errors += 1
	if q < Q*0.7:
		print("ERROR Q %.4f << %.4f %s" % (q, Q, fn))
		errors += 1
	if tc < TC*0.7:
		print("ERROR T %.4f << %.4f %s" % (tc, TC, fn))
		errors += 1

qscore("../test_logs/qscore_BB_seqs_BB11001.log", Q=1, TC=1)
qscore("../test_logs/qscore_BB_seqs_BB11002.log", Q=0.611, TC=0)
qscore("../test_logs/qscore_BB_seqs_BB11004.log", Q=0.674, TC=0.5)
qscore("../test_logs/qscore_BB_seqs_BB11005.log", Q=0.568, TC=0.17)
qscore("../test_logs/qscore_BB_seqs_BB11006.log", Q=0.561, TC=0.375)
qscore("../test_logs/qscore_BB_seqs_BB11007.log", Q=0.837, TC=0.662)
qscore("../test_logs/qscore_BB_seqs_BB11009.log", Q=0.758, TC=0.7)
qscore("../test_logs/qscore_BB_structs_BB11001.log", Q=0.985, TC=0.982)
qscore("../test_logs/qscore_BB_structs_BB11002.log", Q=0.838, TC=0.353)
qscore("../test_logs/qscore_BB_structs_BB11004.log", Q=0.812, TC=0.685)
qscore("../test_logs/qscore_BB_structs_BB11005.log", Q=0.754, TC=0.455)
qscore("../test_logs/qscore_BB_structs_BB11006.log", Q=0.683, TC=0.5)
qscore("../test_logs/qscore_BB_structs_BB11007.log", Q=0.925, TC=0.768)
qscore("../test_logs/qscore_BB_structs_BB11009.log", Q=0.764, TC=0.617)
qscore("../test_logs/qscore_rdrp.log", Q=0.589, TC=0.155)

print("check_results %d errors" % errors)
exit(1 if errors > 0 else 0)
