#include "muscle.h"
#include "tree.h"

#define TRACE	0

// Return false when done
bool PhyEnumEdges(const Tree &tree, PhyEnumEdgeState &ES)
	{
	unsigned uNode1 = UINT_MAX;

	if (!ES.m_bInit)
		{
		if (tree.GetNodeCount() <= 1)
			{
			ES.m_uNodeIndex1 = NULL_NEIGHBOR;
			ES.m_uNodeIndex2 = NULL_NEIGHBOR;
			return false;
			}
		uNode1 = tree.FirstDepthFirstNode();
		ES.m_bInit = true;
		}
	else
		{
		uNode1 = tree.NextDepthFirstNode(ES.m_uNodeIndex1);
		if (NULL_NEIGHBOR == uNode1)
			return false;
		if (tree.IsRooted() && tree.IsRoot(uNode1))
			{
			uNode1 = tree.NextDepthFirstNode(uNode1);
			if (NULL_NEIGHBOR == uNode1)
				return false;
			}
		}
	unsigned uNode2 = tree.GetParent(uNode1);

	ES.m_uNodeIndex1 = uNode1;
	ES.m_uNodeIndex2 = uNode2;
	return true;
	}

bool PhyEnumEdgesR(const Tree &tree, PhyEnumEdgeState &ES)
	{
	unsigned uNode1 = UINT_MAX;

	if (!ES.m_bInit)
		{
		if (tree.GetNodeCount() <= 1)
			{
			ES.m_uNodeIndex1 = NULL_NEIGHBOR;
			ES.m_uNodeIndex2 = NULL_NEIGHBOR;
			return false;
			}
		uNode1 = tree.FirstDepthFirstNodeR();
		ES.m_bInit = true;
		}
	else
		{
		uNode1 = tree.NextDepthFirstNodeR(ES.m_uNodeIndex1);
		if (NULL_NEIGHBOR == uNode1)
			return false;
		if (tree.IsRooted() && tree.IsRoot(uNode1))
			{
			uNode1 = tree.NextDepthFirstNode(uNode1);
			if (NULL_NEIGHBOR == uNode1)
				return false;
			}
		}
	unsigned uNode2 = tree.GetParent(uNode1);

	ES.m_uNodeIndex1 = uNode1;
	ES.m_uNodeIndex2 = uNode2;
	return true;
	}

static void GetLeavesSubtree(const Tree &tree, unsigned uNodeIndex1,
  const unsigned uNodeIndex2, unsigned Leaves[], unsigned *ptruCount)
	{
	if (tree.IsLeaf(uNodeIndex1))
		{
		Leaves[*ptruCount] = uNodeIndex1;
		++(*ptruCount);
		return;
		}

	const unsigned uLeft = tree.GetFirstNeighbor(uNodeIndex1, uNodeIndex2);
	const unsigned uRight = tree.GetSecondNeighbor(uNodeIndex1, uNodeIndex2);
	if (NULL_NEIGHBOR != uLeft)
		GetLeavesSubtree(tree, uLeft, uNodeIndex1, Leaves, ptruCount);
	if (NULL_NEIGHBOR != uRight)
		GetLeavesSubtree(tree, uRight, uNodeIndex1, Leaves, ptruCount);
	}

static void PhyGetLeaves(const Tree &tree, unsigned uNodeIndex1, unsigned uNodeIndex2,
  unsigned Leaves[], unsigned *ptruCount)
	{
	*ptruCount = 0;
	GetLeavesSubtree(tree, uNodeIndex1, uNodeIndex2, Leaves, ptruCount);
	}

bool PhyEnumBiParts(const Tree &tree, PhyEnumEdgeState &ES,
  unsigned Leaves1[], unsigned *ptruCount1,
  unsigned Leaves2[], unsigned *ptruCount2)
	{
	bool bOk = PhyEnumEdges(tree, ES);
	if (!bOk)
		{
		*ptruCount1 = 0;
		*ptruCount2 = 0;
		return false;
		}

// Special case: in a rooted tree, both edges from the root
// give the same bipartition, so skip one of them.
	if (tree.IsRooted() && tree.IsRoot(ES.m_uNodeIndex2)
	  && tree.GetRight(ES.m_uNodeIndex2) == ES.m_uNodeIndex1)
		{
		bOk = PhyEnumEdges(tree, ES);
		if (!bOk)
			return false;
		}

	PhyGetLeaves(tree, ES.m_uNodeIndex1, ES.m_uNodeIndex2, Leaves1, ptruCount1);
	PhyGetLeaves(tree, ES.m_uNodeIndex2, ES.m_uNodeIndex1, Leaves2, ptruCount2);

	if (*ptruCount1 + *ptruCount2 != tree.GetLeafCount())
		Die("PhyEnumBiParts %u + %u != %u",
		  *ptruCount1, *ptruCount2, tree.GetLeafCount());
#if	DEBUG
	{
	for (unsigned i = 0; i < *ptruCount1; ++i)
		{
		if (!tree.IsLeaf(Leaves1[i]))
			Die("PhyEnumByParts: not leaf");
		for (unsigned j = 0; j < *ptruCount2; ++j)
			{
			if (!tree.IsLeaf(Leaves2[j]))
				Die("PhyEnumByParts: not leaf");
			if (Leaves1[i] == Leaves2[j])
				Die("PhyEnumByParts: dupe");
			}
		}
	}
#endif

	return true;
	}

#if	0
void TestBiPart()
	{
	SetListFileName("c:\\tmp\\lobster.log", false);
	Tree tree;
	TextFile fileIn("c:\\tmp\\test.phy");
	tree.FromFile(fileIn);
	tree.LogMe();

	const unsigned uNodeCount = tree.GetNodeCount();
	unsigned *Leaves1 = new unsigned[uNodeCount];
	unsigned *Leaves2 = new unsigned[uNodeCount];

	PhyEnumEdgeState ES;
	bool bDone = false;
	for (;;)
		{
		unsigned uCount1 = UINT_MAX;
		unsigned uCount2 = UINT_MAX;
		bool bOk = PhyEnumBiParts(tree, ES, Leaves1, &uCount1, Leaves2, &uCount2);
		Log("PEBP=%d ES.Init=%d ES.ni1=%d ES.ni2=%d\n",
		  bOk,
		  ES.m_bInit,
		  ES.m_uNodeIndex1,
		  ES.m_uNodeIndex2);
		if (!bOk)
			break;
		Log("\n");
		Log("Part1: ");
		for (unsigned n = 0; n < uCount1; ++n)
			Log(" %d(%s)", Leaves1[n], tree.GetLeafName(Leaves1[n]));
		Log("\n");
		Log("Part2: ");
		for (unsigned n = 0; n < uCount2; ++n)
			Log(" %d(%s)", Leaves2[n], tree.GetLeafName(Leaves2[n]));
		Log("\n");
		}
	}
#endif

static void GetLeavesSubtreeExcluding(const Tree &tree, unsigned uNodeIndex,
  unsigned uExclude, unsigned Leaves[], unsigned *ptruCount)
	{
	if (uNodeIndex == uExclude)
		return;

	if (tree.IsLeaf(uNodeIndex))
		{
		Leaves[*ptruCount] = uNodeIndex;
		++(*ptruCount);
		return;
		}

	const unsigned uLeft = tree.GetLeft(uNodeIndex);
	const unsigned uRight = tree.GetRight(uNodeIndex);
	if (NULL_NEIGHBOR != uLeft)
		GetLeavesSubtreeExcluding(tree, uLeft, uExclude, Leaves, ptruCount);
	if (NULL_NEIGHBOR != uRight)
		GetLeavesSubtreeExcluding(tree, uRight, uExclude, Leaves, ptruCount);
	}

void GetLeavesExcluding(const Tree &tree, unsigned uNodeIndex,
  unsigned uExclude, unsigned Leaves[], unsigned *ptruCount)
	{
	*ptruCount = 0;
	GetLeavesSubtreeExcluding(tree, uNodeIndex, uExclude, Leaves, ptruCount);
	}

void GetInternalNodesInHeightOrder(const Tree &tree, unsigned NodeIndexes[])
	{
	const unsigned uNodeCount = tree.GetNodeCount();
	if (uNodeCount < 3)
		Die("GetInternalNodesInHeightOrder: %u nodes, none are internal",
		  uNodeCount);
	const unsigned uInternalNodeCount = (uNodeCount - 1)/2;
	double *Heights = new double[uInternalNodeCount];

	unsigned uIndex = 0;
	for (unsigned uNodeIndex = 0; uNodeIndex < uNodeCount; ++uNodeIndex)
		{
		if (tree.IsLeaf(uNodeIndex))
			continue;
		NodeIndexes[uIndex] = uNodeIndex;
		Heights[uIndex] = tree.GetNodeHeight(uNodeIndex);
		++uIndex;
		}
	if (uIndex != uInternalNodeCount)
		Die("Internal error: GetInternalNodesInHeightOrder");

// Simple but slow bubble sort (probably don't care about speed here)
	bool bDone = false;
	while (!bDone)
		{
		bDone = true;
		for (unsigned i = 0; i < uInternalNodeCount - 1; ++i)
			{
			if (Heights[i] > Heights[i+1])
				{
				double dTmp = Heights[i];
				Heights[i] = Heights[i+1];
				Heights[i+1] = dTmp;

				unsigned uTmp = NodeIndexes[i];
				NodeIndexes[i] = NodeIndexes[i+1];
				NodeIndexes[i+1] = uTmp;
				bDone = false;
				}
			}
		}
#if	TRACE
	Log("Internal node index     Height\n");
	Log("-------------------   --------\n");
	//    1234567890123456789  123456789
	for (unsigned n = 0; n < uInternalNodeCount; ++n)
		Log("%19u  %9.3f\n", NodeIndexes[n], Heights[n]);
#endif
	delete[] Heights;
	}

void ApplyMinEdgeLength(Tree &tree, double dMinEdgeLength)
	{
	const unsigned uNodeCount = tree.GetNodeCount();
	for (unsigned uNodeIndex = 0; uNodeIndex < uNodeCount; ++uNodeIndex)
		{
		const unsigned uNeighborCount = tree.GetNeighborCount(uNodeIndex);
		for (unsigned n = 0; n < uNeighborCount; ++n)
			{
			const unsigned uNeighborNodeIndex = tree.GetNeighbor(uNodeIndex, n);
			if (!tree.HasEdgeLength(uNodeIndex, uNeighborNodeIndex))
				continue;
			if (tree.GetEdgeLength(uNodeIndex, uNeighborNodeIndex) < dMinEdgeLength)
				tree.SetEdgeLength(uNodeIndex, uNeighborNodeIndex, dMinEdgeLength);
			}
		}
	}
