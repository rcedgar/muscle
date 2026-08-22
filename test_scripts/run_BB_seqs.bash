#!/bin/bash -e

outdir=../test_output/BB_seqs
logdir=../test_logs
mkdir -p $outdir $logdir

for acc in `cat ../test_data/info/BB.accs`
do
	../bin/muscle \
	  -align ../test_data/fa/$acc \
	  -output $outdir/$acc \
	  -log $logdir/BB_seqs.$acc.log
done
