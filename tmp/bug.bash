#!/bin/bash -e

../github_releases/muscle-linux-x86.v5.3 \
	-align uniprotkb_p450_1A1_AND_reviewed_true_2025_09_08.fasta \
	-output p450s.msa.afa \
	-threads 4 \
	-log muscle.log
