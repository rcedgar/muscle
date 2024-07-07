#include "myutils.h"
#include "muscle.h"
#include "kmerdist33.h"

uint KmerDist33::SeqToKmer(const byte *Seq) const
	{
	uint u1 = g_CharToLetterAmino[Seq[0]];
	uint u2 = g_CharToLetterAmino[Seq[1]];
	uint u3 = g_CharToLetterAmino[Seq[2]];

	return u1 + u2*20 + u3*20*20;
	}

void KmerDist33::CountKmers(const byte *Seq, uint L, byte *KmerToCount)
	{
	memset(KmerToCount, 0, DICT_SIZE_33);
	for (uint i = 0; i + 5 < L; ++i)
		{
		uint Kmer = SeqToKmer(Seq+i);
		if (Kmer >= DICT_SIZE_33)
			continue;
		KmerToCount[Kmer] += 1;
		}
	}

uint KmerDist33::GetCommonKmerCount(const byte *KmerToCount1,
  const byte *KmerToCount2) const
	{
	uint Sum = 0;
	for (uint Kmer = 0; Kmer < DICT_SIZE_33; ++Kmer)
		{
		uint n1 = KmerToCount1[Kmer];
		uint n2 = KmerToCount2[Kmer];
		uint MinCount = min(n1, n2);
		Sum += MinCount;
		}
	return Sum;
	}

void KmerDist33::GetDistMx(const MultiSequence &MS,
  vector<vector<float> > &DistMx)
	{
	vector<vector<float> > CommonKmerCountMx;

	const uint SeqCount = MS.GetSeqCount();
	DistMx.resize(SeqCount);
	CommonKmerCountMx.resize(SeqCount);
	for (uint i = 0; i < SeqCount; ++i)
		{
		DistMx[i].resize(SeqCount);
		CommonKmerCountMx[i].resize(SeqCount);
		}

	byte *KmerToCounti = myalloc(byte, DICT_SIZE_33);
	byte *KmerToCountj = myalloc(byte, DICT_SIZE_33);
	uint PairCount = (SeqCount*(SeqCount - 1))/2;
	uint PairIndex = 0;
	for (uint SeqIndexi = 0; SeqIndexi < SeqCount; ++SeqIndexi)
		{
		uint Li;
		const byte *Seqi = MS.GetByteSeq(SeqIndexi, Li);
		CountKmers(Seqi, Li, KmerToCounti);
		float CommonCountii = (float) GetCommonKmerCount(KmerToCounti, KmerToCounti);
		DistMx[SeqIndexi][SeqIndexi] = 0;

		for (uint SeqIndexj = 0; SeqIndexj < SeqIndexi; ++SeqIndexj)
			{
			ProgressStep(PairIndex++, PairCount, "Kmer33 distance");
			uint Lj;
			const byte *Seqj = MS.GetByteSeq(SeqIndexj, Lj);
			CountKmers(Seqj, Lj, KmerToCountj);
			float CommonCountjj = (float) GetCommonKmerCount(KmerToCountj, KmerToCountj);
			float CommonCountij = (float) GetCommonKmerCount(KmerToCounti, KmerToCountj);

			const float d1 = 3.0f*(CommonCountii - CommonCountij)/CommonCountii;
			const float d2 = 3.0f*(CommonCountjj - CommonCountij)/CommonCountjj;
			const float dMin = min(d1, d2);
			DistMx[SeqIndexi][SeqIndexj] = dMin;
			DistMx[SeqIndexj][SeqIndexi] = dMin;
#if 0
			{
			const char *Label1 = MS.GetLabel(SeqIndexi);
			const char *Label2 = MS.GetLabel(SeqIndexj);
			Log("%6.0f  %6.0f  %6.0f  %8.3g  %s, %s\n",
			  CommonCountjj, CommonCountii, CommonCountjj,
			  dMin, Label1, Label2);
			}
#endif
			}
		}

	myfree(KmerToCounti);
	myfree(KmerToCountj);
	}
