#!/bin/bash -e

if [ ! -s ../test_output/rdrp/rdrp_structs.afa ]; then
	echo "SKIP qscore_rdrp: missing ../test_output/rdrp/rdrp_structs.afa"
	exit 0
fi

../bin/muscle \
  -qscore ../test_output/rdrp/rdrp_seqs.afa \
  -ref ../test_output/rdrp/rdrp_structs.afa \
  -bysequence \
  -log ../test_logs/qscore_rdrp.log
