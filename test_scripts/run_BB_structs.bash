#!/bin/bash -e

# Requires STRUCTS fixtures under ../test_data/structs/<acc>/ (or .bca/.files).
# Legacy ../test_data/mega/*.mega is no longer supported.
structs_root=../test_data/structs
if [ ! -d "$structs_root" ]; then
	echo "SKIP run_BB_structs: missing $structs_root (STRUCTS fixtures not checked in)"
	exit 0
fi

outdir=../test_output/BB_structs
logdir=../test_logs
mkdir -p $outdir $logdir

for acc in `cat ../test_data/info/BB.accs`
do
	in="$structs_root/$acc"
	if [ ! -e "$in" ] && [ ! -e "$in.bca" ] && [ ! -e "$in.files" ]; then
		echo "SKIP $acc: no STRUCTS under $structs_root"
		continue
	fi
	[ -e "$in.bca" ] && in="$in.bca"
	[ -e "$in.files" ] && in="$in.files"
	../bin/muscle \
	  -align "$in" \
	  -output $outdir/$acc \
	  -log $logdir/BB_structs.$acc.log
done
