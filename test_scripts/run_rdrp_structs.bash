#!/bin/bash -e

outdir=../test_output/rdrp
logdir=../test_logs
mkdir -p $outdir $logdir ../tmp

cp ../test_data/rdrp.mega.gz ../tmp
gunzip ../tmp/rdrp.mega.gz

../bin/muscle \
  -super7 ../tmp/rdrp/rdrp.mega \
  -guidetreein ../test_data/rdrp/rdrp.newick \
  -output $outdir/rdrp_structs.afa \
  -log ../test_logs/super7_rdrp.log

rm -rf ../tmp
