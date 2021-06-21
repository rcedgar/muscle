#include "muscle.h"
#include "msa.h"
#include "tree.h"
#include "profile.h"
#include "pwpath.h"

#define TRACE	0

// Progressive alignment according to a diffs tree.

void DeleteProgNode(ProgNode &Node)
	{
	delete[] Node.m_Prof;
	delete[] Node.m_EstringL;
	delete[] Node.m_EstringR;

	Node.m_Prof = 0;
	Node.m_EstringL = 0;
	Node.m_EstringR = 0;
	}

static void MakeNode(const MSA &msaIn, const Tree &Diffs, unsigned uDiffsNodeIndex,
   const unsigned IdToDiffsTreeNodeIndex[], ProgNode &Node)
	{
	const unsigned uSeqCount = msaIn.GetSeqCount();

	unsigned *Ids = new unsigned[uSeqCount];

	unsigned uSeqsInDiffCount = 0;
	for (unsigned uId = 0; uId < uSeqCount; ++uId)
		{
		if (IdToDiffsTreeNodeIndex[uId] == uDiffsNodeIndex)
			{
			Ids[uSeqsInDiffCount] = uId;
			++uSeqsInDiffCount;
			}
		}
	if (0 == uSeqsInDiffCount)
		Quit("MakeNode: no seqs in diff");

	MSASubsetByIds(msaIn, Ids, uSeqsInDiffCount, Node.m_MSA);

#if	DEBUG
	ValidateMuscleIds(Node.m_MSA);
#endif

	DeleteGappedCols(Node.m_MSA);
	delete[] Ids;
	}

void RealignDiffs(const MSA &msaIn, const Tree &Diffs,
  const unsigned IdToDiffsTreeNodeIndex[], MSA &msaOut)
	{
	assert(Diffs.IsRooted());

#if	TRACE
	Log("RealignDiffs\n");
	Log("Diff tree:\n");
	Diffs.LogMe();
#endif

	const unsigned uNodeCount = Diffs.GetNodeCount();
	if (uNodeCount%2 == 0)
		Quit("RealignDiffs: Expected odd number of nodes");

	const unsigned uMergeCount = (uNodeCount - 1)/2;

	ProgNode *ProgNodes = new ProgNode[uNodeCount];

	unsigned uJoin = 0;
	//SetProgressDesc("Refine tree");
	uint Counter = 0;
	for (unsigned uDiffsNodeIndex = Diffs.FirstDepthFirstNode();
	  NULL_NEIGHBOR != uDiffsNodeIndex;
	  uDiffsNodeIndex = Diffs.NextDepthFirstNode(uDiffsNodeIndex))
		{
		ProgressStep(Counter++, uNodeCount, "Refine by tree");
		if (Diffs.IsLeaf(uDiffsNodeIndex))
			{
			assert(uDiffsNodeIndex < uNodeCount);
			if (uDiffsNodeIndex >= uNodeCount)
				Quit("TreeNodeIndex=%u NodeCount=%u\n", uDiffsNodeIndex, uNodeCount);

			ProgNode &Node = ProgNodes[uDiffsNodeIndex];
			MakeNode(msaIn, Diffs, uDiffsNodeIndex, IdToDiffsTreeNodeIndex, Node);

			Node.m_uLength = Node.m_MSA.GetColCount();
			}
		else
			{
//			Progress(uJoin, uMergeCount);
			++uJoin;
			const unsigned uMergeNodeIndex = uDiffsNodeIndex;
			ProgNode &Parent = ProgNodes[uMergeNodeIndex];

			const unsigned uLeft = Diffs.GetLeft(uDiffsNodeIndex);
			const unsigned uRight = Diffs.GetRight(uDiffsNodeIndex);

			ProgNode &Node1 = ProgNodes[uLeft];
			ProgNode &Node2 = ProgNodes[uRight];

			PWPath Path;
			AlignTwoMSAs(Node1.m_MSA, Node2.m_MSA, Parent.m_MSA, Path);

#if	TRACE
			{
			Log("Combined:\n");
			Parent.m_MSA.LogMe();
			}
#endif

			Node1.m_MSA.Clear();
			Node2.m_MSA.Clear();
			}
		}
//	ProgressStepsDone();

	unsigned uRootNodeIndex = Diffs.GetRootNodeIndex();
	const ProgNode &RootProgNode = ProgNodes[uRootNodeIndex];
	msaOut.Copy(RootProgNode.m_MSA);

#if	DEBUG
	AssertMSAEqIgnoreCaseAndGaps(msaIn, msaOut);
#endif

	delete[] ProgNodes;
	ProgNodes = 0;
	}
