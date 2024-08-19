#!/bin/bash -e

outdir=../test_output/rdrp
logdir=../test_logs
mkdir -p $outdir $logdir

../bin/muscle \
  -super5 ../test_data/rdrp/rdrp.fa \
  -output $outdir/rdrp_seqs.afa \
  -log ../test_logs/super5_rdrp.log
