#include "myutils.h"
#include "probcons.h"

/////////////////////////////////////////////////////////////////
// DoIterativeRefinement()
//
// Performs a single round of randomized partionining iterative
// refinement.
/////////////////////////////////////////////////////////////////

void DoIterativeRefinement(
  const vector<vector<SparseMatrix*> >& sparseMatrices,
  MultiSequence*& alignment)
	{
	set<int> groupOne, groupTwo;

	// create two separate groups
	for (int i = 0; i < alignment->GetNumSequences(); i++)
		{
		if (rand() % 2)
			groupOne.insert(i);
		else
			groupTwo.insert(i);
		}

	if (groupOne.empty() || groupTwo.empty())
		return;

// project into the two groups
	MultiSequence* groupOneSeqs = alignment->Project(groupOne);
	assert(groupOneSeqs);

	MultiSequence* groupTwoSeqs = alignment->Project(groupTwo);
	assert(groupTwoSeqs);
	delete alignment;

// realign
	alignment = AlignAlignments(groupOneSeqs, groupTwoSeqs, sparseMatrices);

	delete groupOneSeqs;
	delete groupTwoSeqs;
	}
