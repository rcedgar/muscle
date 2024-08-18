#!/bin/bash -e

outdir=../test_output/rdrp
logdir=../test_logs
mkdir -p $outdir $logdir

../bin/muscle \
  -super7 ../test_data/rdrp/rdrp.mega \
  -guidetreein ../test_data/rdrp/rdrp.newick \
  -output $outdir/rdrp_structs.afa
