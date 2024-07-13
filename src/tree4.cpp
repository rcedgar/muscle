#include "muscle.h"
#include "tree.h"
#include <stdio.h>

#define	TRACE	0

void ClusterByHeight(const Tree &tree, double dMaxHeight, unsigned Subtrees[],
  unsigned *ptruSubtreeCount)
	{
	if (!tree.IsRooted())
		Die("ClusterByHeight: requires rooted tree");

#if	TRACE
	Log("ClusterByHeight, max height=%g\n", dMaxHeight);
#endif

	unsigned uSubtreeCount = 0;
	const unsigned uNodeCount = tree.GetNodeCount();
	for (unsigned uNodeIndex = 0; uNodeIndex < uNodeCount; ++uNodeIndex)
		{
		if (tree.IsRoot(uNodeIndex))
			continue;
		unsigned uParent = tree.GetParent(uNodeIndex);
		double dHeight = tree.GetNodeHeight(uNodeIndex);
		double dParentHeight = tree.GetNodeHeight(uParent);

#if	TRACE
		Log("Node %3u  Height %5.2f  ParentHeight %5.2f\n",
		  uNodeIndex, dHeight, dParentHeight);
#endif
		if (dParentHeight > dMaxHeight && dHeight <= dMaxHeight)
			{
			Subtrees[uSubtreeCount] = uNodeIndex;
#if	TRACE
			Log("Subtree[%u]=%u\n", uSubtreeCount, uNodeIndex);
#endif
			++uSubtreeCount;
			}
		}
	*ptruSubtreeCount = uSubtreeCount;
	}

static void ClusterBySubfamCount_Iteration(const Tree &tree, unsigned Subfams[],
  unsigned uCount)
	{
// Find highest child node of current set of subfamilies.
	double dHighestHeight = -1e20;
	int iParentSubscript = -1;

	for (int n = 0; n < (int) uCount; ++n)
		{
		const unsigned uNodeIndex = Subfams[n];
		if (tree.IsLeaf(uNodeIndex))
			continue;

		const unsigned uLeft = tree.GetLeft(uNodeIndex);
		const double dHeightLeft = tree.GetNodeHeight(uLeft);
		if (dHeightLeft > dHighestHeight)
			{
			dHighestHeight = dHeightLeft;
			iParentSubscript = n;
			}

		const unsigned uRight = tree.GetRight(uNodeIndex);
		const double dHeightRight = tree.GetNodeHeight(uRight);
		if (dHeightRight > dHighestHeight)
			{
			dHighestHeight = dHeightRight;
			iParentSubscript = n;
			}
		}

	if (-1 == iParentSubscript)
		Die("CBSFCIter: failed to find highest child");

	const unsigned uNodeIndex = Subfams[iParentSubscript];
	const unsigned uLeft = tree.GetLeft(uNodeIndex);
	const unsigned uRight = tree.GetRight(uNodeIndex);

// Delete parent by replacing with left child
	Subfams[iParentSubscript] = uLeft;

// Append right child to list
	Subfams[uCount] = uRight;

#if	TRACE
	{
	Log("Iter %3u:", uCount);
	for (unsigned n = 0; n < uCount; ++n)
		Log(" %u", Subfams[n]);
	Log("\n");
	}
#endif
	}

// Divide a tree containing N leaves into k families by
// cutting the tree at a horizontal line at some height.
// Each internal node defines a height for the cut, 
// considering all internal nodes enumerates all distinct
// cuts. Visit internal nodes in decreasing order of height.
// Visiting the node corresponds to moving the horizontal
// line down to cut the tree at the height of that node.
// We consider the cut to be "infinitestimally below"
// the node, so the effect is to remove the current node 
// from the list of subfamilies and add its two children.
// We must visit a parent before its children (so care may
// be needed to handle zero edge lengths properly).
// We assume that N is small, and write dumb O(N^2) code.
// More efficient strategies are possible for large N
// by maintaining a list of nodes sorted by height.
void ClusterBySubfamCount(const Tree &tree, unsigned uSubfamCount,
  unsigned Subfams[], unsigned *ptruSubfamCount)
	{
	const unsigned uNodeCount = tree.GetNodeCount();
	const unsigned uLeafCount = (uNodeCount + 1)/2;

// Special case: empty tree
	if (0 == uNodeCount)
		{
		*ptruSubfamCount = 0;
		return;
		}

// Special case: more subfamilies than leaves
	if (uSubfamCount >= uLeafCount)
		{
		for (unsigned n = 0; n < uLeafCount; ++n)
			Subfams[n] = n;
		*ptruSubfamCount = uLeafCount;
		return;
		}

// Initialize list of subfamilies to be root
	Subfams[0] = tree.GetRootNodeIndex();

// Iterate
	for (unsigned i = 1; i < uSubfamCount; ++i)
		ClusterBySubfamCount_Iteration(tree, Subfams, i);
	
	*ptruSubfamCount = uSubfamCount;
	}

static void GetLeavesRecurse(const Tree &tree, unsigned uNodeIndex,
  unsigned Leaves[], unsigned &uLeafCount /* in-out */)
	{
	if (tree.IsLeaf(uNodeIndex))
		{
		Leaves[uLeafCount] = uNodeIndex;
		++uLeafCount;
		return;
		}

	const unsigned uLeft = tree.GetLeft(uNodeIndex);
	const unsigned uRight = tree.GetRight(uNodeIndex);

	GetLeavesRecurse(tree, uLeft, Leaves, uLeafCount);
	GetLeavesRecurse(tree, uRight, Leaves, uLeafCount);
	}

void GetLeaves(const Tree &tree, unsigned uNodeIndex, unsigned Leaves[],
  unsigned *ptruLeafCount)
	{
	unsigned uLeafCount = 0;
	GetLeavesRecurse(tree, uNodeIndex, Leaves, uLeafCount);
	*ptruLeafCount = uLeafCount;
	}

void Tree::PruneTree(const Tree &tree, unsigned Subfams[],
  unsigned uSubfamCount, const char *LabelPrefix, vector<string> &Labels)
	{
	if (!tree.IsRooted())
		Die("Tree::PruneTree: requires rooted tree");
	Labels.clear();

	Clear();

	m_uNodeCount = 2*uSubfamCount - 1;
	InitCache(m_uNodeCount);

	const unsigned uUnprunedNodeCount = tree.GetNodeCount();

	unsigned *uUnprunedToPrunedIndex = new unsigned[uUnprunedNodeCount];
	unsigned *uPrunedToUnprunedIndex = new unsigned[m_uNodeCount];

	for (unsigned n = 0; n < uUnprunedNodeCount; ++n)
		uUnprunedToPrunedIndex[n] = NULL_NEIGHBOR;

	for (unsigned n = 0; n < m_uNodeCount; ++n)
		uPrunedToUnprunedIndex[n] = NULL_NEIGHBOR;

// Create mapping between unpruned and pruned node indexes
	unsigned uInternalNodeIndex = uSubfamCount;
	for (unsigned uSubfamIndex = 0; uSubfamIndex < uSubfamCount; ++uSubfamIndex)
		{
		unsigned uUnprunedNodeIndex = Subfams[uSubfamIndex];
		uUnprunedToPrunedIndex[uUnprunedNodeIndex] = uSubfamIndex;
		uPrunedToUnprunedIndex[uSubfamIndex] = uUnprunedNodeIndex;
		for (;;)
			{
			uUnprunedNodeIndex = tree.GetParent(uUnprunedNodeIndex);
			if (tree.IsRoot(uUnprunedNodeIndex))
				break;

		// Already visited this node?
			if (NULL_NEIGHBOR != uUnprunedToPrunedIndex[uUnprunedNodeIndex])
				break;

			uUnprunedToPrunedIndex[uUnprunedNodeIndex] = uInternalNodeIndex;
			uPrunedToUnprunedIndex[uInternalNodeIndex] = uUnprunedNodeIndex;

			++uInternalNodeIndex;
			}
		}

	const unsigned uUnprunedRootIndex = tree.GetRootNodeIndex();
	uUnprunedToPrunedIndex[uUnprunedRootIndex] = uInternalNodeIndex;
	uPrunedToUnprunedIndex[uInternalNodeIndex] = uUnprunedRootIndex;

#if	TRACE
	{
	Log("Pruned to unpruned:\n");
	for (unsigned i = 0; i < m_uNodeCount; ++i)
		Log(" [%u]=%u", i, uPrunedToUnprunedIndex[i]);
	Log("\n");
	Log("Unpruned to pruned:\n");
	for (unsigned i = 0; i < uUnprunedNodeCount; ++i)
		{
		unsigned n = uUnprunedToPrunedIndex[i];
		if (n != NULL_NEIGHBOR)
			Log(" [%u]=%u", i, n);
		}
	Log("\n");
	}
#endif

	if (uInternalNodeIndex != m_uNodeCount - 1)
		Die("Tree::PruneTree, Internal error");

// Nodes 0, 1 ... are the leaves
	for (unsigned uSubfamIndex = 0; uSubfamIndex < uSubfamCount; ++uSubfamIndex)
		{
		string Label;
		Ps(Label, "%s%u", LabelPrefix, uSubfamIndex);
		m_ptrName[uSubfamIndex] = mystrsave(Label.c_str());
		Labels.push_back(Label);
		}

	for (unsigned uPrunedNodeIndex = uSubfamCount; uPrunedNodeIndex < m_uNodeCount;
	  ++uPrunedNodeIndex)
		{
		unsigned uUnprunedNodeIndex = uPrunedToUnprunedIndex[uPrunedNodeIndex];

		const unsigned uUnprunedLeft = tree.GetLeft(uUnprunedNodeIndex);
		const unsigned uUnprunedRight = tree.GetRight(uUnprunedNodeIndex);

		const unsigned uPrunedLeft = uUnprunedToPrunedIndex[uUnprunedLeft];
		const unsigned uPrunedRight = uUnprunedToPrunedIndex[uUnprunedRight];

		const double dLeftLength =
		  tree.GetEdgeLength(uUnprunedNodeIndex, uUnprunedLeft);
		const double dRightLength =
		  tree.GetEdgeLength(uUnprunedNodeIndex, uUnprunedRight);

		m_uNeighbor2[uPrunedNodeIndex] = uPrunedLeft;
		m_uNeighbor3[uPrunedNodeIndex] = uPrunedRight;

		m_dEdgeLength1[uPrunedLeft] = dLeftLength;
		m_dEdgeLength1[uPrunedRight] = dRightLength;

		m_uNeighbor1[uPrunedLeft] = uPrunedNodeIndex;
		m_uNeighbor1[uPrunedRight] = uPrunedNodeIndex;

		m_bHasEdgeLength1[uPrunedLeft] = true;
		m_bHasEdgeLength1[uPrunedRight] = true;

		m_dEdgeLength2[uPrunedNodeIndex] = dLeftLength;
		m_dEdgeLength3[uPrunedNodeIndex] = dRightLength;

		m_bHasEdgeLength2[uPrunedNodeIndex] = true;
		m_bHasEdgeLength3[uPrunedNodeIndex] = true;
		}

	m_uRootNodeIndex = uUnprunedToPrunedIndex[uUnprunedRootIndex];

	m_bRooted = true;

	Validate();

	delete[] uUnprunedToPrunedIndex;
	}

void LeafIndexesToIds(const Tree &tree, const unsigned Leaves[], unsigned uCount,
  unsigned Ids[])
	{
	for (unsigned n = 0; n < uCount; ++n)
		Ids[n] = tree.GetLeafId(Leaves[n]);
	}
