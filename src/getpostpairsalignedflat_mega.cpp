#include "muscle.h"
#include "locallock.h"

float GetPostPairsAlignedFlat_mega(const string &aProgressStr,
  const MultiSequence &MSA1, const MultiSequence &MSA2,
  const vector<uint> &SeqIndexes1, const vector<uint> &SeqIndexes2, 
  vector<MySparseMx *> &SparsePosts)
	{
	string ProgressStr = aProgressStr;
	if (SIZE(ProgressStr) > 20)
		ProgressStr = ProgressStr.substr(0, 20);

	const uint SeqCount1 = MSA1.GetSeqCount();
	const uint SeqCount2 = MSA2.GetSeqCount();
	const uint PairCount = SIZE(SeqIndexes1);
	asserta(SIZE(SeqIndexes2) == PairCount);
	asserta(SparsePosts.empty());

// Allocate here to avoid race condition with push_back() in loop
	SparsePosts.resize(PairCount);

	int PairCounter = 0;
	uint ThreadCount = GetRequestedThreadCount();
	float SumEA = 0;
#pragma omp parallel for num_threads(ThreadCount)
	for (int PairIndex = 0; PairIndex < (int) PairCount; ++PairIndex)
		{
		uint Min = min(SeqCount1, SeqCount2);
		uint Max = max(SeqCount1, SeqCount2);
		Lock();
		ProgressStep(PairCounter++, PairCount,
		  "%s [%u x %u, %u pairs]",
		  ProgressStr.c_str(), Min, Max, PairCount);
		Unlock();

		uint SeqIndex1 = SeqIndexes1[PairIndex];
		uint SeqIndex2 = SeqIndexes2[PairIndex];
		asserta(SeqIndex1 < SeqCount1);
		asserta(SeqIndex2 < SeqCount2);

		const Sequence *gapped_seq1 = MSA1.GetSequence(SeqIndex1);
		const Sequence *gapped_seq2 = MSA2.GetSequence(SeqIndex2);
		Sequence *seq1 = gapped_seq1->CopyDeleteGaps();
		Sequence *seq2 = gapped_seq2->CopyDeleteGaps();
		const byte *ByteSeq1 = seq1->GetBytePtr();
		const byte *ByteSeq2 = seq2->GetBytePtr();
		const uint L1 = seq1->GetLength();
		const uint L2 = seq2->GetLength();

		float *Fwd = AllocFB(L1, L2);
		float *Bwd = AllocFB(L1, L2);
		float *Post = AllocPost(L1, L2);

		CalcFwdFlat(ByteSeq1, L1, ByteSeq2, L2, Fwd);
		CalcBwdFlat(ByteSeq1, L1, ByteSeq2, L2, Bwd);

		DeleteSequence(seq1);
		DeleteSequence(seq2);

		CalcPostFlat(Fwd, Bwd, L1, L2, Post);
		delete Fwd;
		delete Bwd;

		float *DPRows = AllocDPRows(L1, L2);
		char *TB = AllocTB(L1, L2);

		string Path;
		float Score = CalcAlnFlat(Post, L1, L2, DPRows, TB, Path);
		delete DPRows;
		delete TB;

		MySparseMx *SparsePost = new MySparseMx;
		asserta(SparsePost);
		SparsePost->FromPost(Post, L1, L2);
		SparsePosts[PairIndex] = SparsePost;
		delete Post;

		float EA = Score/min(L1, L2);
		Lock();
		SumEA += EA;
		Unlock();
		}
	float AvgEA = SumEA/PairCount;
	return AvgEA;
	}
