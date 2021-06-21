#include "myutils.h"
#include "muscle.h"
#include "probcons.h"

void ProgressLogInputSummary(const string &FileName, const MultiSequence &Seqs);

void CalcEADistMx(FILE *f, MultiSequence* sequences,
  vector<vector<float> > &DistMx)
	{
	DistMx.clear();
	const uint SeqCount = sequences->GetSeqCount();
	DistMx.resize(SeqCount);
	for (uint i = 0; i < SeqCount; ++i)
		{
		DistMx[i].resize(SeqCount, 0);
		DistMx[i][i] = 1;
		}

	uint PairIndex = 0;
	vector<pair<uint, uint> > Pairs;
	for (uint SeqIndex1 = 0; SeqIndex1 < SeqCount; SeqIndex1++)
		for (uint SeqIndex2 = SeqIndex1 + 1; SeqIndex2 < SeqCount; SeqIndex2++)
			Pairs.push_back(pair<uint, uint>(SeqIndex1, SeqIndex2));
	uint PairCount2 = SIZE(Pairs);
	uint PairCount = (SeqCount * (SeqCount - 1)) / 2;
	asserta(PairCount == PairCount2);

// all-vs-all pairwise alignments for posterior probability matrices
	unsigned ThreadCount = GetRequestedThreadCount();
//	ProgressLog("%u fwd-bwd threads\n", ThreadCount);
	omp_lock_t Lock;
	omp_init_lock(&Lock);
	uint PairCounter = 0;
	float SumEA = 0;
	float MinEA = 1;
	float MaxEA = 0;
#pragma omp parallel for num_threads(ThreadCount)
	for (int PairIndex = 0; PairIndex < (int) PairCount; ++PairIndex)
		{
		const pair<uint, uint>& Pair = Pairs[PairIndex];
		int SeqIndex1 = Pair.first;
		int SeqIndex2 = Pair.second;

		Sequence* seq1 = sequences->GetSequence(SeqIndex1);
		Sequence* seq2 = sequences->GetSequence(SeqIndex2);

		const char *Label1 = seq1->label.c_str();
		const char *Label2 = seq2->label.c_str();

		int L1 = seq1->GetLength();
		int L2 = seq2->GetLength();

		omp_set_lock(&Lock);
		double MeanEA = (PairCounter == 0 ? 0 : SumEA/PairCounter);
		ProgressStep(PairCounter++, PairCount,
		  "%u consensus seqs, EE min %.2g mean %.2g max %.2g", 
		  SeqCount, 1 - MinEA, 1 - MeanEA, 1 - MaxEA);
		omp_unset_lock(&Lock);

		vector<float>* forward =
		  PairHMM::ComputeForwardMatrix(seq1, seq2); assert(forward);

		vector<float>* backward =
		  PairHMM::ComputeBackwardMatrix(seq1, seq2); assert(backward);

		vector<float>* posterior =
		  PairHMM::ComputePosteriorMatrix(seq1, seq2, *forward, *backward);

		pair<vector<char>*, float> alignment =
		  PairHMM::ComputeAlignment(L1, L2, *posterior);

		float EA = alignment.second / min(L1, L2);
		omp_set_lock(&Lock);
		DistMx[SeqIndex1][SeqIndex2] = EA;
		DistMx[SeqIndex2][SeqIndex1] = EA;
		if (f != 0)
			fprintf(f, "%s\t%s\t%.4g\n", Label1, Label2, EA);
		SumEA += EA;
		MinEA = min(MinEA, EA);
		MaxEA = max(MaxEA, EA);
		omp_unset_lock(&Lock);

		delete alignment.first;
		delete posterior;
		delete forward;
		delete backward;
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
