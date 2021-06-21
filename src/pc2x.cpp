#include "myutils.h"
#include "probcons.h"

void cmd_pc2x()
	{
	const string &MSAFileName1 = opt(pc2x);
	const string &MSAFileName2 = opt(input2);

	int TargetPairCount = 1024;
	if (optset_paircount)
		TargetPairCount = int(opt(paircount));

	//SetProbconsParams();
	//const PairHMM HMM(initDistrib, gapOpen, gapExtend,
	//  emitPairs, emitSingle);
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

	int AllPairCount = SeqCount1*SeqCount2;
	int PairCount = UINT_MAX;
	vector<int> SeqIndexes1;
	vector<int> SeqIndexes2;
	if (AllPairCount < TargetPairCount*3/2)
		{
		for (int i = 0; i < SeqCount1; ++i)
			for (int j = 0; j < SeqCount2; ++j)
				{
				SeqIndexes1.push_back(i);
				SeqIndexes2.push_back(j);
				}
		asserta(SIZE(SeqIndexes1) == AllPairCount);
		asserta(SIZE(SeqIndexes2) == AllPairCount);
		PairCount = AllPairCount;
		}
	else
		{
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
		PairCount = SIZE(PairSet);
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
	ProgressLog("%u pairs (target %u)\n", PairCount, TargetPairCount);
//	vector<float> distances(PairCount);
	vector<SparseMatrix*> sparseMatrices(PairCount);

	omp_lock_t Lock;
	omp_init_lock(&Lock);
	int PairCounter = 0;
	uint ThreadCount = GetRequestedThreadCount();
	double SumDistance = 0;
	ProgressLog("%u fwd-bwd threads\n", ThreadCount);
#pragma omp parallel for num_threads(ThreadCount)
	for (int PairIndex = 0; PairIndex < (int) PairCount; ++PairIndex)
		{
		omp_set_lock(&Lock);
		ProgressStep(PairCounter++, PairCount, "Pairwise fwd-bwd");
		omp_unset_lock(&Lock);

		int SeqIndex1 = SeqIndexes1[PairIndex];
		int SeqIndex2 = SeqIndexes2[PairIndex];
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
		sparseMatrices[PairIndex] =
		  new SparseMatrix(SeqLength1, SeqLength2, *posterior);

	// perform the pairwise sequence alignment
		pair<vector<char>*, float> alignment = PairHMM::ComputeAlignment(
		  SeqLength1, SeqLength2, *posterior);

	// compute "expected accuracy" distance
		float distance = alignment.second / min(seq1->GetLength(), seq2->GetLength());
		omp_set_lock(&Lock);
		//distances[PairIndex] = distance;
		SumDistance += distance;
		omp_unset_lock(&Lock);

		delete seq1;
		delete seq2;
		delete alignment.first;
		delete posterior;
		delete forward;
		delete backward;
		}
	
	double MeanDistance = SumDistance/PairCount;
	ProgressLog("MeanDistance=%.4g\n", MeanDistance);
	
	vector<float>* posterior = PairHMM::BuildPosterior3(&MSA1, &MSA2, 
	  SeqIndexes1, SeqIndexes2, sparseMatrices);

	pair<vector<char>*, float> alignment = 
	  PairHMM::ComputeAlignment(ColCount1, ColCount2, *posterior);
	const vector<char> &Path = *alignment.first;
	float Score = alignment.second;

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

	if (optset_output)
		result->WriteMFA(opt(output));

	delete result;
	delete alignment.first;
	}
