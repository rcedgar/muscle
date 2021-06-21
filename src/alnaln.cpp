#include "myutils.h"
#include "probcons.h"

MultiSequence* AlignAlignments(MultiSequence* align1, MultiSequence* align2,
  const vector<vector<SparseMatrix*> > &SparseMatrices)
	{
	vector<float>* posterior = 
	  PairHMM::BuildPosterior(align1, align2, SparseMatrices);
	pair<vector<char>*, float> alignment;

	alignment = PairHMM::ComputeAlignment(align1->GetSequence(0)->GetLength(),
	  align2->GetSequence(0)->GetLength(), *posterior);
	float Score = alignment.second;
	delete posterior;

	MultiSequence* result = new MultiSequence();
	for (int i = 0; i < align1->GetNumSequences(); i++)
		result->AddSequence(align1->GetSequence(i)->AddGaps(alignment.first, 'X'));
	for (int i = 0; i < align2->GetNumSequences(); i++)
		result->AddSequence(align2->GetSequence(i)->AddGaps(alignment.first, 'Y'));

	delete alignment.first;
	return result;
	}
