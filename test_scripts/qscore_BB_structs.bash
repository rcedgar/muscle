#!/bin/bash -e

outdir=../test_output/BB_structs
if [ ! -d "$outdir" ] || [ -z "$(ls -A "$outdir" 2>/dev/null)" ]; then
	echo "SKIP qscore_BB_structs: no BB_structs outputs"
	exit 0
fi

logdir=../test_logs
mkdir -p $outdir $logdir

for acc in `cat ../test_data/info/BB.accs`
do
	if [ ! -s "$outdir/$acc" ]; then
		echo "SKIP qscore $acc: missing $outdir/$acc"
		continue
	fi
	../bin/muscle \
	  -qscore ../test_output/BB_structs/$acc \
	  -ref ../test_data/ref_alns/$acc \
	  -bysequence \
	  -log ../test_logs/qscore_BB_structs_$acc.log
done
