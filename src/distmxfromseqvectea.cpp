#include "myutils.h"
#include "muscle.h"
#include "seqvect.h"
#include "probcons.h"

Sequence *SeqToSequence(const Seq &inseq, uint Index);

MultiSequence *SeqVectToMultiSequence(const SeqVect &SV)
	{
	MultiSequence *MS = new MultiSequence;
	asserta(MS != 0);
	const uint SeqCount = SV.GetSeqCount();
	for (uint i = 0; i < SeqCount; ++i)
		{
		const Seq &seq1 = SV.GetSeq(i);
		Sequence *seq2 = SeqToSequence(seq1, i);
		MS->AddSequence(seq2);
		}
	return MS;
	}

// TODO@@ Lots of shared code with RunMPC()
void DistMxFromSeqVect_EA(const SeqVect &SV, vector<vector<float> > &DistMx,
  vector<string> &Labels)
	{
	DistMx.clear();
	Labels.clear();

	MultiSequence *sequences = SeqVectToMultiSequence(SV);

	DistMx.clear();
	Labels.clear();

	const uint SeqCount = SV.GetSeqCount();
	DistMx.resize(SeqCount);
	for (uint i = 0; i < SeqCount; ++i)
		{
		DistMx[i].resize(SeqCount, FLT_MAX);

		string Label;
		SV.GetLabel(i, Label);
		Labels.push_back(Label);
		}

	uint PairIndex = 0;
	vector<pair<uint, uint> > Pairs;
	for (uint SeqIndex1 = 0; SeqIndex1 < SeqCount; SeqIndex1++)
		for (uint SeqIndex2 = SeqIndex1 + 1; SeqIndex2 < SeqCount; SeqIndex2++)
			Pairs.push_back(pair<uint, uint>(SeqIndex1, SeqIndex2));

	uint PairCount2 = SIZE(Pairs);
	uint PairCount = (SeqCount * (SeqCount - 1)) / 2;
	asserta(PairCount == PairCount2);

	vector<vector<SparseMatrix*> > sparseMatrices(SeqCount, vector<SparseMatrix*>(SeqCount, NULL));

// all-vs-all pairwise alignments for posterior probability matrices
	unsigned ThreadCount = GetRequestedThreadCount();
	ProgressLog("distmx/EA %u fwd-bwd threads\n", ThreadCount);
	omp_lock_t Lock;
	omp_init_lock(&Lock);
	uint PairCounter = 0;
#pragma omp parallel for num_threads(ThreadCount)
	for (int PairIndex = 0; PairIndex < (int) PairCount; ++PairIndex)
		{
		const pair<uint, uint>& Pair = Pairs[PairIndex];
		int SeqIndex1 = Pair.first;
		int SeqIndex2 = Pair.second;
		Sequence* seq1 = sequences->GetSequence(SeqIndex1);
		Sequence* seq2 = sequences->GetSequence(SeqIndex2);
		omp_set_lock(&Lock);
		ProgressStep(PairCounter++, PairCount, "distmx/EA fwd-bwd");
		omp_unset_lock(&Lock);

	// compute forward and backward probabilities
		vector<float>* forward = PairHMM::ComputeForwardMatrix(seq1, seq2); assert(forward);
		vector<float>* backward = PairHMM::ComputeBackwardMatrix(seq1, seq2); assert(backward);

	// compute posterior probability matrix
		vector<float>* posterior = PairHMM::ComputePosteriorMatrix(seq1, seq2, *forward, *backward); assert(posterior);

	// compute sparse representations
		sparseMatrices[SeqIndex1][SeqIndex2] =
		  new SparseMatrix(seq1->GetLength(), seq2->GetLength(), *posterior);
		sparseMatrices[SeqIndex2][SeqIndex1] = 0;

		pair<vector<char>*, float> alignment = PairHMM::ComputeAlignment(
		  seq1->GetLength(),
		  seq2->GetLength(),
		  *posterior);

	// "expected accuracy" distance
		float EA = alignment.second / min(seq1->GetLength(), seq2->GetLength());
		DistMx[SeqIndex1][SeqIndex2] = EA;
		DistMx[SeqIndex2][SeqIndex1] = EA;

		delete alignment.first;
		delete posterior;
		delete forward;
		delete backward;
		}

  // delete sparse matrices
	PairCounter = 0;
	for (uint SeqIndex1 = 0; SeqIndex1 < SeqCount - 1; SeqIndex1++)
		{
		for (uint SeqIndex2 = SeqIndex1 + 1; SeqIndex2 < SeqCount; SeqIndex2++)
			{
			ProgressStep(PairCounter++, PairCount, "distmx/EA Delete sparse matrices");
			SparseMatrix *SM12 = sparseMatrices[SeqIndex1][SeqIndex2];
			SparseMatrix *SM21 = sparseMatrices[SeqIndex2][SeqIndex1];
			if (SM12 != 0)
				delete SM12;
			if (SM21 != SM12 && SM21 != 0)
				delete SM21;
			}
		}
	sparseMatrices.clear();
	}
