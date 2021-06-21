#include "muscle.h"
#include <math.h>
#include "tree.h"
#include "seqvect.h"
#include "profile.h"
#include "msa.h"
#include "pwpath.h"
#include "distfunc.h"

#define TRACE 0

void ProgressiveAlign(const SeqVect &v, const Tree &GuideTree, MSA &a)
	{
	assert(GuideTree.IsRooted());

#if	TRACE
	Log("GuideTree:\n");
	GuideTree.LogMe();
#endif

	const unsigned uSeqCount = v.Length();
	const unsigned uNodeCount = 2*uSeqCount - 1;

	ProgNode *ProgNodes = new ProgNode[uNodeCount];

	unsigned uJoin = 0;
	unsigned uTreeNodeIndex = GuideTree.FirstDepthFirstNode();
	//SetProgressDesc("Align node");
	do
		{
		if (GuideTree.IsLeaf(uTreeNodeIndex))
			{
			if (uTreeNodeIndex >= uNodeCount)
				Quit("TreeNodeIndex=%u NodeCount=%u\n", uTreeNodeIndex, uNodeCount);
			ProgNode &Node = ProgNodes[uTreeNodeIndex];
			unsigned uId = GuideTree.GetLeafId(uTreeNodeIndex);
			if (uId >= uSeqCount)
				Quit("Seq index out of range");
			const Seq &s = *(v[uId]);
			Node.m_MSA.FromSeq(s);
			Node.m_MSA.SetSeqId(0, uId);
			Node.m_uLength = Node.m_MSA.GetColCount();
			}
		else
			{
			//Progress(uJoin, uSeqCount - 1);
			ProgressStep(uJoin, uSeqCount - 1, "Progressive alignment");
			++uJoin;

			const unsigned uMergeNodeIndex = uTreeNodeIndex;
			ProgNode &Parent = ProgNodes[uMergeNodeIndex];

			const unsigned uLeft = GuideTree.GetLeft(uTreeNodeIndex);
			const unsigned uRight = GuideTree.GetRight(uTreeNodeIndex);

			ProgNode &Node1 = ProgNodes[uLeft];
			ProgNode &Node2 = ProgNodes[uRight];

			PWPath Path;
			AlignTwoMSAs(Node1.m_MSA, Node2.m_MSA, Parent.m_MSA, Path);
			Parent.m_uLength = Parent.m_MSA.GetColCount();

			Node1.m_MSA.Clear();
			Node2.m_MSA.Clear();
			}
		uTreeNodeIndex = GuideTree.NextDepthFirstNode(uTreeNodeIndex);
		}
	while (NULL_NEIGHBOR != uTreeNodeIndex);
	ProgressStepsDone();

	unsigned uRootNodeIndex = GuideTree.GetRootNodeIndex();
	const ProgNode &RootProgNode = ProgNodes[uRootNodeIndex];
	a.Copy(RootProgNode.m_MSA);

	delete[] ProgNodes;
	ProgNodes = 0;
	}
