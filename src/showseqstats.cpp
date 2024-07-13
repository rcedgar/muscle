#include "muscle.h"

void ShowSeqStats(const MultiSequence &InputSeqs)
	{
	const uint InputSeqCount = InputSeqs.GetSeqCount();
	double MeanSeqLength = InputSeqs.GetMeanSeqLength();
	uint MaxSeqLength = InputSeqs.GetMaxSeqLength();
	uint MinSeqLength = InputSeqs.GetMinSeqLength();
	ProgressLog("Input: %u seqs, avg length %.0f, max %u, min %u\n\n",
	  InputSeqCount, MeanSeqLength, MaxSeqLength, MinSeqLength);
	if (MaxSeqLength > 100000)
		Die("Too long, not appropriate for global alignment");
	if (MaxSeqLength > 20000)
		Warning("Long sequences, likey to crash / not globally alignable, max length is ~21k");
	}
