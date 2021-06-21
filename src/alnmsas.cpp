#include "myutils.h"
#include "probcons.h"
#include "sort.h"

#define TRACE	0
#define TOPSORT	0

static void GetAllPairs(uint SeqCount1, uint SeqCount2,
  vector<int> &SeqIndexes1, vector<int> &SeqIndexes2)
	{
	SeqIndexes1.clear();
	SeqIndexes2.clear();
	for (int i = 0; i < (int) SeqCount1; ++i)
		{
		for (int j = 0; j < (int) SeqCount2; ++j)
			{
			SeqIndexes1.push_back(i);
			SeqIndexes2.push_back(j);
			}
		}
	}

static void GetPairs(uint SeqCount1, uint SeqCount2, uint TargetPairCount,
  vector<int> &SeqIndexes1, vector<int> &SeqIndexes2)
	{
	SeqIndexes1.clear();
	SeqIndexes2.clear();

	uint AllPairCount = SeqCount1*SeqCount2;
	if (TargetPairCount == UINT_MAX || AllPairCount < TargetPairCount*3/2)
		{
		GetAllPairs(SeqCount1, SeqCount2, SeqIndexes1, SeqIndexes2);
		return;
		}

	set<pair<int, int> > PairSet;
	const int MaxCounter = TargetPairCount*10;
	int Counter = 0;
	while (Counter++ < MaxCounter && (int) SIZE(PairSet) < TargetPairCount)
		{
		int i = randu32()%SeqCount1;
		int j = randu32()%SeqCount2;
		if (i == j)
			continue;
		pair<int, int> Pair(i, j);
		PairSet.insert(Pair);
		}

	uint PairCount = SIZE(PairSet);
	asserta(PairCount > TargetPairCount/2);
	for (set<pair<int, int> >::const_iterator p = PairSet.begin();
		p != PairSet.end(); ++p)
		{
		int SeqIndex1 = p->first;
		int SeqIndex2 = p->second;
		SeqIndexes1.push_back(SeqIndex1);
		SeqIndexes2.push_back(SeqIndex2);
		}
	}

void _AssertSeqsEq(const char *FileName, uint LineNr,
  const MultiSequence &MSA1, const MultiSequence &MSA2)
	{
	const uint SeqCount1 = MSA1.GetSeqCount();
	for (uint SeqIndex1 = 0; SeqIndex1 < SeqCount1; ++SeqIndex1)
		{
		const Sequence *Seq1 = MSA1.GetSequence((int) SeqIndex1);
		const string &Label = Seq1->label;
		uint SeqIndex2 = MSA2.GetSeqIndex(Label);
		const Sequence *Seq2 = MSA2.GetSequence((int) SeqIndex2);

		uint GSI1 = Seq1->GetGSI();
		uint GSI2 = Seq2->GetGSI();

		Sequence *uSeq1 = Seq1->DeleteGaps();
		Sequence *uSeq2 = Seq2->DeleteGaps();
		int Length1 = uSeq1->GetLength();
		int Length2 = uSeq2->GetLength();

		const vector<char> &v1 = *(uSeq1->data);
		const vector<char> &v2 = *(uSeq2->data);
		if (v1 != v2 || GSI1 != GSI2)
			{
			Log("\n");
			Log("AssertSeqsEq >%s\n", Label.c_str());
			Log("GI1 %u, GI2 %u\n", GSI1, GSI2);
			Log("Seq1[%d]  ", Length1);
			for (int i = 1; i < Length1; ++i)
				Log("%c", v1[i]);
			Log("\n");
			Log("Seq2[%d]  ", Length2);
			for (int i = 1; i < Length2; ++i)
				Log("%c", v2[i]);
			Log("\n");
			Die("AssertSeqsEq %s:%u", FileName, LineNr);
			}

		delete uSeq1;
		delete uSeq2;
		}
	}

float AlignMSAs(const string &aProgressStr,
  const MultiSequence &MSA1, 
  const MultiSequence &MSA2,
  uint TargetPairCount,
  vector<char> &Path)
	{
	string ProgressStr = aProgressStr;
	if (SIZE(ProgressStr) > 20)
		ProgressStr = ProgressStr.substr(0, 20);
	const int SeqCount1 = MSA1.GetNumSequences();
	const int SeqCount2 = MSA2.GetNumSequences();
	asserta(SeqCount1 > 0);
	asserta(SeqCount2 > 0);

	asserta(MSA1.IsAligned());
	asserta(MSA2.IsAligned());

	const int ColCount1 = MSA1.GetSequence(0)->GetLength();
	const int ColCount2 = MSA2.GetSequence(0)->GetLength();

	vector<int> SeqIndexes1;
	vector<int> SeqIndexes2;
	GetPairs(SeqCount1, SeqCount2, TargetPairCount,
	  SeqIndexes1, SeqIndexes2);
	const uint PairCount = SIZE(SeqIndexes1);
	asserta(SIZE(SeqIndexes2) == PairCount);
	Log("%s: %u pairs (target %u)\n",
	  ProgressStr.c_str(), PairCount, TargetPairCount);

	vector<SparseMatrix*> sparseMatrices(PairCount);

	omp_lock_t Lock;
	omp_init_lock(&Lock);
	int PairCounter = 0;
	uint ThreadCount = GetRequestedThreadCount();
	float SumExpectedAccuracy = 0;
#if TOPSORT
	vector<float> EAs;
#endif
#pragma omp parallel for num_threads(ThreadCount)
	for (int PairIndex = 0; PairIndex < (int) PairCount; ++PairIndex)
		{
		omp_set_lock(&Lock);
		ProgressStep(PairCounter++, PairCount,
		  "%s fwd-bwd %u threads %u pairs", ProgressStr.c_str(), ThreadCount, PairCount);
		omp_unset_lock(&Lock);

		int SeqIndex1 = SeqIndexes1[PairIndex];
		int SeqIndex2 = SeqIndexes2[PairIndex];
		asserta(SeqIndex1 < SeqCount1);
		asserta(SeqIndex2 < SeqCount2);

		const Sequence *gapped_seq1 = MSA1.GetSequence(SeqIndex1);
		const Sequence *gapped_seq2 = MSA2.GetSequence(SeqIndex2);
		Sequence *seq1 = gapped_seq1->DeleteGaps();
		Sequence *seq2 = gapped_seq2->DeleteGaps();

	// compute forward and backward probabilities
		vector<float>* forward = PairHMM::ComputeForwardMatrix(seq1, seq2);
		vector<float>* backward = PairHMM::ComputeBackwardMatrix(seq1, seq2);
		asserta(forward != 0);
		assert(backward != 0);

		vector<float>* posterior = 
		  PairHMM::ComputePosteriorMatrix(seq1, seq2, *forward, *backward);
		asserta(posterior != 0);

		int SeqLength1 = seq1->GetLength();
		int SeqLength2 = seq2->GetLength();
		sparseMatrices[PairIndex] =
		  new SparseMatrix(SeqLength1, SeqLength2, *posterior);

		pair<vector<char>*, float> alignment =
		  PairHMM::ComputeAlignment(SeqLength1, SeqLength2, *posterior);
		float Score = alignment.second;

		int L1 = seq1->GetLength();
		int L2 = seq2->GetLength();
		float ExpectedAccuracy = Score/min(L1, L2);
#if TOPSORT
		EAs.push_back(ExpectedAccuracy);
#endif
		omp_set_lock(&Lock);
		SumExpectedAccuracy += ExpectedAccuracy;
		omp_unset_lock(&Lock);

		delete seq1;
		delete seq2;
		delete alignment.first;
		delete posterior;
		delete forward;
		delete backward;
		}
	float MeanExpectedAccuracy = SumExpectedAccuracy/PairCount;

#if TOPSORT
	{
	vector<uint> Order(PairCount);
	QuickSortOrderDesc(EAs.data(), PairCount, Order.data());
	vector<SparseMatrix*> TopSparseMatrices;
	vector<int> TopSeqIndexes1;
	vector<int> TopSeqIndexes2;
	for (uint i = 0; i < PairCount/2; ++i)
		{
		uint k = Order[i];
		int SeqIndex1 = SeqIndexes1[k];
		int SeqIndex2 = SeqIndexes2[k];
		TopSeqIndexes1.push_back(SeqIndex1);
		TopSeqIndexes2.push_back(SeqIndex2);
		TopSparseMatrices.push_back(sparseMatrices[k]);
		}
	}
	vector<float>* posterior = PairHMM::BuildPosterior3(&MSA1, &MSA2,
	  TopSeqIndexes1, TopSeqIndexes2, TopSparseMatrices);
#else	
	vector<float>* posterior = PairHMM::BuildPosterior3(&MSA1, &MSA2,
	  SeqIndexes1, SeqIndexes2, sparseMatrices);
#endif

	for (uint i = 0; i < PairCount; ++i)
		delete sparseMatrices[i];
	sparseMatrices.clear();

	pair<vector<char>*, float> alignment = 
	  PairHMM::ComputeAlignment(ColCount1, ColCount2, *posterior);
	Path = *alignment.first;
	float AlnScore = alignment.second;
	const uint PathLength = SIZE(Path);

	delete posterior;
	delete alignment.first;
	return MeanExpectedAccuracy;
	}
