#include "myutils.h"
#include "muscle.h"
#include "locallock.h"

void ProgressLogInputSummary(const string &FileName, const MultiSequence &Seqs);

void CalcEADistMx(FILE *f, MultiSequence* sequences,
  vector<vector<float> > &DistMx, vector<MySparseMx * > *SparsePostVec)
	{
	DistMx.clear();
	const uint SeqCount = sequences->GetSeqCount();
	DistMx.resize(SeqCount);
	for (uint i = 0; i < SeqCount; ++i)
		{
		DistMx[i].resize(SeqCount, 0);
		DistMx[i][i] = 1;
		}
	if (SparsePostVec != 0)
		asserta(SIZE(*SparsePostVec) == 0);

	vector<uint> SeqIndexes1;
	vector<uint> SeqIndexes2;
	GetAllPairs(SeqCount, SeqIndexes1, SeqIndexes2);
	uint PairCount = SIZE(SeqIndexes1);
	asserta(SIZE(SeqIndexes1) == PairCount);
	uint PairCount2 = (SeqCount * (SeqCount - 1)) / 2;
	asserta(PairCount == PairCount2);

// all-vs-all pairwise alignments for posterior probability matrices
	unsigned ThreadCount = GetRequestedThreadCount();
	uint PairCounter = 0;
	float SumEA = 0;
#pragma omp parallel for num_threads(ThreadCount)
	for (int PairIndex = 0; PairIndex < (int) PairCount; ++PairIndex)
		{
		uint SeqIndex1 = SeqIndexes1[PairIndex];
		uint SeqIndex2 = SeqIndexes2[PairIndex];

		const Sequence* seq1 = sequences->GetSequence(SeqIndex1);
		const Sequence* seq2 = sequences->GetSequence(SeqIndex2);

		const string &Label1 = seq1->m_Label;
		const string &Label2 = seq2->m_Label;

		Lock();
		double MeanEA = (PairCounter == 0 ? 0 : SumEA/PairCounter);
		ProgressStep(PairCounter++, PairCount,
		  "%u consensus seqs, mean EE %.2g", SeqCount, 1 - MeanEA);
		Unlock();

		string Path;
		float EA;
		if (SparsePostVec == 0)
			EA = AlignPairFlat(Label1, Label2, Path);
		else
			{
			MySparseMx *SparsePost = new MySparseMx;
			EA = AlignPairFlat_SparsePost(Label1, Label2, Path, SparsePost);
			SparsePostVec->push_back(SparsePost);
			}

		Lock();
		DistMx[SeqIndex1][SeqIndex2] = EA;
		DistMx[SeqIndex2][SeqIndex1] = EA;
		if (f != 0)
			fprintf(f, "%s\t%s\t%.4g\n", Label1.c_str(), Label2.c_str(), EA);
		SumEA += EA;
		Unlock();
		}
	}

void cmd_eadistmx()
	{
	const string &InputFileName = opt(eadistmx);
	asserta(optset_output);
	FILE *f = CreateStdioFile(opt(output));

	MultiSequence* sequences = new MultiSequence();
	assert(sequences);
	sequences->LoadMFA(InputFileName, true);
	ProgressLogInputSummary(InputFileName, *sequences);

	InitProbcons();

	vector<vector<float> > DistMx;
	CalcEADistMx(f, sequences, DistMx);
	CloseStdioFile(f);
	}
