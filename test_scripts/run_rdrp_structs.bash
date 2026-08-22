#!/bin/bash -e

# Requires STRUCTS (e.g. ../test_data/rdrp/rdrp.bca). Text .mega is no longer supported.
structs=../test_data/rdrp/rdrp.bca
if [ ! -s "$structs" ]; then
	echo "SKIP run_rdrp_structs: missing $structs"
	exit 0
fi

outdir=../test_output/rdrp
logdir=../test_logs
mkdir -p $outdir $logdir

../bin/muscle \
  -super7 "$structs" \
  -guidetreein ../test_data/rdrp/rdrp.newick \
  -output $outdir/rdrp_structs.afa \
  -log ../test_logs/super7_rdrp.log
