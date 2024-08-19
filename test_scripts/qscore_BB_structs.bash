#!/bin/bash -e

outdir=../test_output/BB_structs
logdir=../test_logs
mkdir -p $outdir $logdir

for acc in `cat ../test_data/info/BB.accs`
do
	../bin/muscle \
	  -qscore ../test_output/BB_structs/$acc \
	  -ref ../test_data/ref_alns/$acc \
	  -bysequence \
	  -log ../test_logs/qscore_BB_structs_$acc.log
done
