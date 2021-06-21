#include "muscle.h"
#include "tree.h"
#include "msa.h"

/***
Compute weights by the CLUSTALW method.
Thompson, Higgins and Gibson (1994), CABIOS (10) 19-29;
see also CLUSTALW paper.

Weights are computed from the edge lengths of a rooted tree.

Define the strength of an edge to be its length divided by the number
of leaves under that edge. The weight of a sequence is then the sum
of edge strengths on the path from the root to the leaf.

Example.

        0.2
       -----A     0.1
	 -x         ------- B     0.7
	   --------y           ----------- C
	    0.3     ----------z
                    0.4    -------------- D
                                 0.8

Edge	Length	Leaves	Strength
----	-----	------	--------
xy		0.3		3		0.1
xA		0.2		1		0.2
yz		0.4		2		0.2
yB		0.1		1		0.1
zC		0.7		1		0.7
zD		0.8		1		0.8

Leaf	Path		Strengths			Weight
----	----		---------			------
A		xA			0.2					0.2
B		xy-yB		0.1 + 0.1			0.2
C		xy-yz-zC	0.1 + 0.2 + 0.7		1.0
D		xy-yz-zD	0.1 + 0.2 + 0.8		1.1

***/

#define TRACE 0

static unsigned CountLeaves(const Tree &tree, unsigned uNodeIndex,
  unsigned LeavesUnderNode[])
	{
	if (tree.IsLeaf(uNodeIndex))
		{
		LeavesUnderNode[uNodeIndex] = 1;
		return 1;
		}

	const unsigned uLeft = tree.GetLeft(uNodeIndex);
	const unsigned uRight = tree.GetRight(uNodeIndex);
	const unsigned uRightCount = CountLeaves(tree, uRight, LeavesUnderNode);
	const unsigned uLeftCount = CountLeaves(tree, uLeft, LeavesUnderNode);
	const unsigned uCount = uRightCount + uLeftCount;
	LeavesUnderNode[uNodeIndex] = uCount;
	return uCount;
	}

void CalcClustalWWeights(const Tree &tree, WEIGHT Weights[])
	{
#if	TRACE
	Log("CalcClustalWWeights\n");
	tree.LogMe();
#endif

	const unsigned uLeafCount = tree.GetLeafCount();
	if (0 == uLeafCount)
		return;
	else if (1 == uLeafCount)
		{
		Weights[0] = (WEIGHT) 1.0;
		return;
		}
	else if (2 == uLeafCount)
		{
		Weights[0] = (WEIGHT) 0.5;
		Weights[1] = (WEIGHT) 0.5;
		return;
		}

	if (!tree.IsRooted())
		Quit("CalcClustalWWeights requires rooted tree");

	const unsigned uNodeCount = tree.GetNodeCount();
	unsigned *LeavesUnderNode = new unsigned[uNodeCount];
	memset(LeavesUnderNode, 0, uNodeCount*sizeof(unsigned));

	const unsigned uRootNodeIndex = tree.GetRootNodeIndex();
	unsigned uLeavesUnderRoot = CountLeaves(tree, uRootNodeIndex, LeavesUnderNode);
	if (uLeavesUnderRoot != uLeafCount)
		Quit("WeightsFromTreee: Internal error, root count %u %u",
		  uLeavesUnderRoot, uLeafCount);

#if	TRACE
	Log("Node  Leaves    Length  Strength\n");
	Log("----  ------  --------  --------\n");
	//    1234  123456  12345678  12345678
#endif

	double *Strengths = new double[uNodeCount];
	for (unsigned uNodeIndex = 0; uNodeIndex < uNodeCount; ++uNodeIndex)
		{
		if (tree.IsRoot(uNodeIndex))
			{
			Strengths[uNodeIndex] = 0.0;
			continue;
			}
		const unsigned uParent = tree.GetParent(uNodeIndex);
		const double dLength = tree.GetEdgeLength(uNodeIndex, uParent);
		const unsigned uLeaves = LeavesUnderNode[uNodeIndex];
		const double dStrength = dLength / (double) uLeaves;
		Strengths[uNodeIndex] = dStrength;
#if	TRACE
		Log("%4u  %6u  %8g  %8g\n", uNodeIndex, uLeaves, dLength, dStrength);
#endif
		}

#if	TRACE
	Log("\n");
	Log("                 Seq  Path..Weight\n");
	Log("--------------------  ------------\n");
#endif
	for (unsigned n = 0; n < uLeafCount; ++n)
		{
		const unsigned uLeafNodeIndex = tree.LeafIndexToNodeIndex(n);
#if	TRACE
		Log("%20.20s  %4u ", tree.GetLeafName(uLeafNodeIndex), uLeafNodeIndex);
#endif
		if (!tree.IsLeaf(uLeafNodeIndex))
			Quit("CalcClustalWWeights: leaf");

		double dWeight = 0;
		unsigned uNode = uLeafNodeIndex;
		while (!tree.IsRoot(uNode))
			{
			dWeight += Strengths[uNode];
			uNode = tree.GetParent(uNode);
#if	TRACE
			Log("->%u(%g)", uNode, Strengths[uNode]);
#endif
			}
		if (dWeight < 0.0001)
			{
#if	TRACE
			Log("zero->one");
#endif
			dWeight = 1.0;
			}
		Weights[n] = (WEIGHT) dWeight;
#if	TRACE
		Log(" = %g\n", dWeight);
#endif
		}

	delete[] Strengths;
	delete[] LeavesUnderNode;

	Normalize(Weights, uLeafCount);
	}

void MSA::SetClustalWWeights(const Tree &tree)
	{
	const unsigned uSeqCount = GetSeqCount();
	const unsigned uLeafCount = tree.GetLeafCount();

	WEIGHT *Weights = new WEIGHT[uSeqCount];

	CalcClustalWWeights(tree, Weights);

	for (unsigned n = 0; n < uLeafCount; ++n)
		{
		const WEIGHT w = Weights[n];
		const unsigned uLeafNodeIndex = tree.LeafIndexToNodeIndex(n);
		const unsigned uId = tree.GetLeafId(uLeafNodeIndex);
		const unsigned uSeqIndex = GetSeqIndex(uId);
#if	DEBUG
		if (GetSeqName(uSeqIndex) != tree.GetLeafName(uLeafNodeIndex))
			Quit("MSA::SetClustalWWeights: names don't match");
#endif
		SetSeqWeight(uSeqIndex, w);
		}
	NormalizeWeights((WEIGHT) 1.0);

	delete[] Weights;
	}
