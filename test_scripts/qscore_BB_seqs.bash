#!/bin/bash -e

outdir=../test_output/BB_seqs
logdir=../test_logs
mkdir -p $outdir $logdir

for acc in `cat ../test_data/info/BB.accs`
do
	../bin//muscle \
	  -qscore ../test_output/BB_seqs/$acc \
	  -ref ../test_data/ref_alns/$acc \
	  -log ../test_logs/qscore_BB_seqs_$acc.log
done
