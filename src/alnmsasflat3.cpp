#include "muscle.h"
#include "pprog.h"

float PProg::AlignMSAsFlat3(const string &ProgressStr,
  const MultiSequence &MSA1, const MultiSequence &MSA2,
  const vector<MySparseMx *> &SparseMxVec,
  uint Index1, uint Index2,
  uint TargetPairCount, string &Path)
	{
	const uint SeqCount1 = MSA1.GetNumSequences();
	const uint SeqCount2 = MSA2.GetNumSequences();
	asserta(SeqCount1 > 0);
	asserta(SeqCount2 > 0);

	asserta(MSA1.IsAligned());
	asserta(MSA2.IsAligned());

	const uint ColCount1 = MSA1.GetColCount();
	const uint ColCount2 = MSA2.GetColCount();

	vector<uint> SeqIndexes1;
	vector<uint> SeqIndexes2;
	GetPairs(SeqCount1, SeqCount2, TargetPairCount,
	  SeqIndexes1, SeqIndexes2);
	const uint PairCount = SIZE(SeqIndexes1);
	asserta(SIZE(SeqIndexes2) == PairCount);

	vector<MySparseMx *> SparseMxs;
	float AvgEA = GetPostPairsAlignedFlat(ProgressStr, MSA1, MSA2,
	  SeqIndexes1, SeqIndexes2, SparseMxs);

	const uint L1 = ColCount1;
	const uint L2 = ColCount2;

	float *Post = AllocPost(L1, L2);
	CalcPosteriorFlat3(MSA1, MSA2, SeqIndexes1, SeqIndexes2, SparseMxs, Post);

	for (uint i = 0; i < PairCount; ++i)
		delete SparseMxs[i];
	SparseMxs.clear();

	float *DPRows = AllocDPRows(L1, L2);
	char *TB = AllocTB(L1, L2);

	CalcAlnFlat(Post, ColCount1, ColCount2, DPRows, TB, Path);

	delete Post;
	delete DPRows;
	delete TB;

	return AvgEA;
	}
