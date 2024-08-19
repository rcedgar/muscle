#!/bin/bash -e

outdir=../test_output/BB_structs
logdir=../test_logs
mkdir -p $outdir $logdir

for acc in `cat ../test_data/info/BB.accs`
do
	../bin/muscle \
	  -align ../test_data/mega/$acc.mega \
	  -output $outdir/$acc \
	  -log $logdir/BB_structs.$acc.log
done
