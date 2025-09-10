#!/bin/bash -e

../github_releases/muscle-linux-x86.v5.3 \
	-threads 4 \
	-log super5.log \
	-super5 uniprotkb_p450_1A1_AND_reviewed_true_2025_09_08.fasta \
	-output p450s.msa.afa
