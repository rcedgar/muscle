#include "myutils.h"
#include "probcons.h"

/////////////////////////////////////////////////////////////////
// ProcessTree()
//
// Process the tree recursively.  Returns the aligned sequences
// corresponding to a node or leaf of the tree.
/////////////////////////////////////////////////////////////////

static uint g_NodeCount;
static uint g_NodeCounter;

static MultiSequence* ProcessTree_Recursive(const TreeNode* tree, MultiSequence* sequences,
  const vector<vector<SparseMatrix*> >& sparseMatrices)
	{
	if (g_NodeCounter < g_NodeCount)
		ProgressStep(g_NodeCounter++, g_NodeCount, "Progressive (%u/%u)",
		  g_NodeCounter, g_NodeCount);

	MultiSequence* result = 0;

	if (tree->GetSequenceLabel() == -1)
		{
		MultiSequence* alignLeft =
		  ProcessTree_Recursive(tree->GetLeftChild(), sequences, sparseMatrices);
		MultiSequence* alignRight =
		  ProcessTree_Recursive(tree->GetRightChild(), sequences, sparseMatrices);

		assert(alignLeft);
		assert(alignRight);

		result = AlignAlignments(alignLeft, alignRight, sparseMatrices);
		assert(result);

		delete alignLeft;
		delete alignRight;
		}
	else
		{
		result = new MultiSequence(); assert(result);
		result->AddSequence(sequences->GetSequence(tree->GetSequenceLabel())->Clone());
		}

	asserta(result != 0);
	return result;
	}

MultiSequence* ProcessTree(const TreeNode* tree, MultiSequence* sequences,
  const vector<vector<SparseMatrix*> >& sparseMatrices)
	{
	const uint SeqCount = sequences->GetSeqCount();
	g_NodeCount = 2*SeqCount - 1;
	g_NodeCounter = 0;
	MultiSequence *result = ProcessTree_Recursive(tree, sequences, sparseMatrices);
	if (g_NodeCounter + 1 < g_NodeCount)
		ProgressStep(g_NodeCount - 1, g_NodeCount, "Progressive");
	Progress("%u progressive nodes completed\n", g_NodeCounter);
	return result;
	}
