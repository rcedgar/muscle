#pragma once

#include "alpha.h"

class MultiSequence;

static const uint DICT_SIZE_33 = 20*20*20;

class KmerDist33
	{
public:
	uint SeqToKmer(const byte *Seq) const;
	void GetDistMx(const MultiSequence &MS,
	  vector<vector<float> > &DistMx);
	void CountKmers(const byte *Seq, uint L,
	  byte *KmerToCount);
	uint GetCommonKmerCount(const byte *KmerToCount1,
	  const byte *KmerToCount2) const;
	};
