#include "myutils.h"
#include "probcons.h"

#define TRACE	0

void cmd_pc2()
	{
	const string &MSAFileName1 = opt(pc2);
	const string &MSAFileName2 = opt(input2);

	InitProbcons();

	MultiSequence MSA1;
	MultiSequence MSA2;
	Progress("Reading %s ...", MSAFileName1.c_str());
	MSA1.LoadMFA(MSAFileName1, false);
	Progress("done.\n");

	Progress("Reading %s ...", MSAFileName2.c_str());
	MSA2.LoadMFA(MSAFileName2, false);
	Progress("done.\n");

	const int SeqCount1 = MSA1.GetNumSequences();
	const int SeqCount2 = MSA2.GetNumSequences();
	asserta(SeqCount1 > 0);
	asserta(SeqCount2 > 0);

	asserta(MSA1.IsAligned());
	asserta(MSA2.IsAligned());

	const int ColCount1 = MSA1.GetSequence(0)->GetLength();
	const int ColCount2 = MSA2.GetSequence(0)->GetLength();

	const int PairCount = SeqCount1*SeqCount2;

	vector<vector<float>> distances(SeqCount1);
	vector<vector<SparseMatrix*> > sparseMatrices(SeqCount1);
	for (int SeqIndex1 = 0; SeqIndex1 < SeqCount1; ++SeqIndex1)
		{
		distances[SeqIndex1].resize(SeqCount2, FLT_MAX);
		sparseMatrices[SeqIndex1].resize(SeqCount2, 0);
		}

	omp_lock_t Lock;
	omp_init_lock(&Lock);
	int PairCounter = 0;
	uint ThreadCount = GetRequestedThreadCount();
	ProgressLog("%u fwd-bwd threads\n", ThreadCount);
#pragma omp parallel for num_threads(ThreadCount)
	for (int PairIndex = 0; PairIndex < (int) PairCount; ++PairIndex)
		{
		omp_set_lock(&Lock);
		ProgressStep(PairCounter++, PairCount, "Pairwise fwd-bwd");
		omp_unset_lock(&Lock);

		int SeqIndex1 = PairIndex/SeqCount2;
		int SeqIndex2  = PairIndex%SeqCount2;
		asserta(SeqIndex1 < SeqCount1);
		asserta(SeqIndex2 < SeqCount2);

		Sequence *gapped_seq1 = MSA1.GetSequence(SeqIndex1);
		Sequence *gapped_seq2 = MSA2.GetSequence(SeqIndex2);
		Sequence *seq1 = gapped_seq1->DeleteGaps();
		Sequence *seq2 = gapped_seq2->DeleteGaps();

	// compute forward and backward probabilities
		vector<float>* forward = PairHMM::ComputeForwardMatrix(seq1, seq2);
		vector<float>* backward = PairHMM::ComputeBackwardMatrix(seq1, seq2);
		asserta(forward != 0);
		assert(backward != 0);

	// compute posterior probability matrix
		vector<float>* posterior = 
		  PairHMM::ComputePosteriorMatrix(seq1, seq2, *forward, *backward);
		asserta(posterior != 0);

	// compute sparse representations
		int SeqLength1 = seq1->GetLength();
		int SeqLength2 = seq2->GetLength();
		sparseMatrices[SeqIndex1][SeqIndex2] =
		  new SparseMatrix(SeqLength1, SeqLength2, *posterior);

	// perform the pairwise sequence alignment
		pair<vector<char>*, float> alignment = PairHMM::ComputeAlignment(
		  SeqLength1, SeqLength2, *posterior);

	// compute "expected accuracy" distance
		float distance = alignment.second / min(seq1->GetLength(), seq2->GetLength());
		distances[SeqIndex1][SeqIndex2] = distance;

		delete seq1;
		delete seq2;
		delete alignment.first;
		delete posterior;
		delete forward;
		delete backward;
		}

#if TRACE
	{
	Log("\n");
	Log("Distances\n");
	for (int SeqIndex1 = 0; SeqIndex1 < SeqCount1; ++SeqIndex1)
		{
		const Sequence *Seq1 = MSA1.GetSequence(SeqIndex1);
		const string &Label1 = Seq1->label;
		const vector<float> &Row = distances[SeqIndex1];
		for (int SeqIndex2 = 0; SeqIndex2 < SeqCount2; ++SeqIndex2)
			{
			const Sequence *Seq2 = MSA2.GetSequence(SeqIndex2);
			const string &Label2 = Seq2->label;
			double d = Row[SeqIndex2];
			Log("%8.8s  %8.8s  %8.3g\n",
			  Label1.c_str(), Label2.c_str(), d);
			}
		}
	}
#endif
	
	vector<float>* posterior =
	  PairHMM::BuildPosterior2(&MSA1, &MSA2, sparseMatrices);

	pair<vector<char>*, float> alignment = 
	  PairHMM::ComputeAlignment(ColCount1, ColCount2, *posterior);
	const vector<char> &Path = *alignment.first;
	float Score = alignment.second;

#if TRACE
	Log("\n");
	Log("Path=");
	for (uint i = 0; i < SIZE(Path); ++i)
		Log("%c", Path[i]);
	Log("\n");
#endif

	delete posterior;
	posterior = 0;

// now build final alignment
	MultiSequence* result = new MultiSequence();
	for (int SeqIndex = 0; SeqIndex < MSA1.GetNumSequences(); ++SeqIndex)
		{
		const Sequence *Seq1 = MSA1.GetSequence(SeqIndex);
		Sequence *AlignedSeq1 = Seq1->AddGaps(&Path, 'X');
		result->AddSequence(AlignedSeq1);
		}

	for (int SeqIndex = 0; SeqIndex < MSA2.GetNumSequences(); ++SeqIndex)
		{
		const Sequence *Seq2 = MSA2.GetSequence(SeqIndex);
		Sequence *AlignedSeq2 = Seq2->AddGaps(&Path, 'Y');
		result->AddSequence(AlignedSeq2);
		}

	result->WriteMFA(opt(output));

	delete result;
	delete alignment.first;
	}
