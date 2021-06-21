#include "muscle.h"
#include "tree.h"
#include "edgelist.h"

#define TRACE	0

struct EdgeInfo
	{
	EdgeInfo()
		{
		m_bSet = false;
		}
// Is data in this structure valid (i.e, has been set)?
	bool m_bSet;

// Node at start of this edge
	unsigned m_uNode1;

// Node at end of this edge
	unsigned m_uNode2;

// Maximum distance from Node2 to a leaf
	double m_dMaxDistToLeaf;

// Sum of distances from Node2 to all leaves under Node2
	double m_dTotalDistToLeaves;

// Next node on path from Node2 to most distant leaf
	unsigned m_uMaxStep;

// Most distant leaf from Node2 (used for debugging only)
	unsigned m_uMostDistantLeaf;

// Number of leaves under Node2
	unsigned m_uLeafCount;
	};

static void RootByMidLongestSpan(const Tree &tree, EdgeInfo **EIs,
  unsigned *ptruNode1, unsigned *ptruNode2,
  double *ptrdLength1, double *ptrdLength2);
static void RootByMinAvgLeafDist(const Tree &tree, EdgeInfo **EIs,
  unsigned *ptruNode1, unsigned *ptruNode2,
  double *ptrdLength1, double *ptrdLength2);

static void ListEIs(EdgeInfo **EIs, unsigned uNodeCount)
	{
	Log("Node1  Node2  MaxDist  TotDist  MostDist  LeafCount  Step\n");
	Log("-----  -----  -------  -------  --------  ---------  ----\n");
	//    12345  12345  1234567  1234567  12345678  123456789

	for (unsigned uNode = 0; uNode < uNodeCount; ++uNode)
		for (unsigned uNeighbor = 0; uNeighbor < 3; ++uNeighbor)
			{
			const EdgeInfo &EI = EIs[uNode][uNeighbor];
			if (!EI.m_bSet)
				continue;
			Log("%5u  %5u  %7.3g  %7.3g  %8u  %9u",
			  EI.m_uNode1,
			  EI.m_uNode2,
			  EI.m_dMaxDistToLeaf,
			  EI.m_dTotalDistToLeaves,
			  EI.m_uMostDistantLeaf,
			  EI.m_uLeafCount);
			if (NULL_NEIGHBOR != EI.m_uMaxStep)
				Log("  %4u", EI.m_uMaxStep);
			Log("\n");
			}
	}

static void CalcInfo(const Tree &tree, unsigned uNode1, unsigned uNode2, EdgeInfo **EIs)
	{
	const unsigned uNeighborIndex = tree.GetNeighborSubscript(uNode1, uNode2);
	EdgeInfo &EI = EIs[uNode1][uNeighborIndex];
	EI.m_uNode1 = uNode1;
	EI.m_uNode2 = uNode2;

	if (tree.IsLeaf(uNode2))
		{
		EI.m_dMaxDistToLeaf = 0;
		EI.m_dTotalDistToLeaves = 0;
		EI.m_uMaxStep = NULL_NEIGHBOR;
		EI.m_uMostDistantLeaf = uNode2;
		EI.m_uLeafCount = 1;
		EI.m_bSet = true;
		return;
		}

	double dMaxDistToLeaf = -1e29;
	double dTotalDistToLeaves = 0.0;
	unsigned uLeafCount = 0;
	unsigned uMostDistantLeaf = NULL_NEIGHBOR;
	unsigned uMaxStep = NULL_NEIGHBOR;

	const unsigned uNeighborCount = tree.GetNeighborCount(uNode2);
	for (unsigned uSub = 0; uSub < uNeighborCount; ++uSub)
		{
		const unsigned uNode3 = tree.GetNeighbor(uNode2, uSub);
		if (uNode3 == uNode1)
			continue;
		const EdgeInfo &EINext = EIs[uNode2][uSub];
		if (!EINext.m_bSet)
			Quit("CalcInfo: internal error, dist %u->%u not known",
				uNode2, uNode3);


		uLeafCount += EINext.m_uLeafCount;

		const double dEdgeLength = tree.GetEdgeLength(uNode2, uNode3);
		const double dTotalDist = EINext.m_dTotalDistToLeaves +
		  EINext.m_uLeafCount*dEdgeLength;
		dTotalDistToLeaves += dTotalDist;

		const double dDist = EINext.m_dMaxDistToLeaf + dEdgeLength;
		if (dDist > dMaxDistToLeaf)
			{
			dMaxDistToLeaf = dDist;
			uMostDistantLeaf = EINext.m_uMostDistantLeaf;
			uMaxStep = uNode3;
			}
		}
	if (NULL_NEIGHBOR == uMaxStep || NULL_NEIGHBOR == uMostDistantLeaf ||
	  0 == uLeafCount)
		Quit("CalcInfo: internal error 2");

	const double dThisDist = tree.GetEdgeLength(uNode1, uNode2);
	EI.m_dMaxDistToLeaf = dMaxDistToLeaf;
	EI.m_dTotalDistToLeaves = dTotalDistToLeaves;
	EI.m_uMaxStep = uMaxStep;
	EI.m_uMostDistantLeaf = uMostDistantLeaf;
	EI.m_uLeafCount = uLeafCount;
	EI.m_bSet = true;
	}

static bool Known(const Tree &tree, EdgeInfo **EIs, unsigned uNodeFrom,
  unsigned uNodeTo)
	{
	const unsigned uSub = tree.GetNeighborSubscript(uNodeFrom, uNodeTo);
	return EIs[uNodeFrom][uSub].m_bSet;
	}

static bool AllKnownOut(const Tree &tree, EdgeInfo **EIs, unsigned uNodeFrom,
  unsigned uNodeTo)
	{
	const unsigned uNeighborCount = tree.GetNeighborCount(uNodeTo);
	for (unsigned uSub = 0; uSub < uNeighborCount; ++uSub)
		{
		unsigned uNeighborIndex = tree.GetNeighbor(uNodeTo, uSub);
		if (uNeighborIndex == uNodeFrom)
			continue;
		if (!EIs[uNodeTo][uSub].m_bSet)
			return false;
		}
	return true;
	}

void FindRoot(const Tree &tree, unsigned *ptruNode1, unsigned *ptruNode2,
  double *ptrdLength1, double *ptrdLength2,
  ROOT RootMethod)
	{
#if	TRACE
	tree.LogMe();
#endif
	if (tree.IsRooted())
		Quit("FindRoot: tree already rooted");

	const unsigned uNodeCount = tree.GetNodeCount();
	const unsigned uLeafCount = tree.GetLeafCount();

	if (uNodeCount < 2)
		Quit("Root: don't support trees with < 2 edges");

	EdgeInfo **EIs = new EdgeInfo *[uNodeCount];
	for (unsigned uNodeIndex = 0; uNodeIndex < uNodeCount; ++uNodeIndex)
		EIs[uNodeIndex] = new EdgeInfo[3];

	EdgeList Edges;
	for (unsigned uNodeIndex = 0; uNodeIndex < uNodeCount; ++uNodeIndex)
		if (tree.IsLeaf(uNodeIndex))
			{
			unsigned uParent = tree.GetNeighbor1(uNodeIndex);
			Edges.Add(uParent, uNodeIndex);
			}

#if	TRACE
	Log("Edges: ");
	Edges.LogMe();
#endif

// Main loop: iterate until all distances known
	double dAllMaxDist = -1e20;
	unsigned uMaxFrom = NULL_NEIGHBOR;
	unsigned uMaxTo = NULL_NEIGHBOR;
	for (;;)
		{
		EdgeList NextEdges;

#if	TRACE
		Log("\nTop of main loop\n");
		Log("Edges: ");
		Edges.LogMe();
		Log("MDs:\n");
		ListEIs(EIs, uNodeCount);
#endif

	// For all edges
		const unsigned uEdgeCount = Edges.GetCount();
		if (0 == uEdgeCount)
			break;
		for (unsigned n = 0; n < uEdgeCount; ++n)
			{
			unsigned uNodeFrom;
			unsigned uNodeTo;
			Edges.GetEdge(n, &uNodeFrom, &uNodeTo);

			CalcInfo(tree, uNodeFrom, uNodeTo, EIs);
#if	TRACE
			Log("Edge %u -> %u\n", uNodeFrom, uNodeTo);
#endif
			const unsigned uNeighborCount = tree.GetNeighborCount(uNodeFrom);
			for (unsigned i = 0; i < uNeighborCount; ++i)
				{
				const unsigned uNeighborIndex = tree.GetNeighbor(uNodeFrom, i);
				if (!Known(tree, EIs, uNeighborIndex, uNodeFrom) &&
				  AllKnownOut(tree, EIs, uNeighborIndex, uNodeFrom))
					NextEdges.Add(uNeighborIndex, uNodeFrom);
				}
			}
		Edges.Copy(NextEdges);
		}

#if	TRACE
	ListEIs(EIs, uNodeCount);
#endif

	switch (RootMethod)
		{
	case ROOT_MidLongestSpan:
		RootByMidLongestSpan(tree, EIs, ptruNode1, ptruNode2,
		  ptrdLength1, ptrdLength2);
		break;

	case ROOT_MinAvgLeafDist:
		RootByMinAvgLeafDist(tree, EIs, ptruNode1, ptruNode2,
		  ptrdLength1, ptrdLength2);
		break;

	default:
		Quit("Invalid RootMethod=%d", RootMethod);
		}

	for (unsigned uNodeIndex = 0; uNodeIndex < uNodeCount; ++uNodeIndex)
		delete[] EIs[uNodeIndex];
	delete[] EIs;
	}

static void RootByMidLongestSpan(const Tree &tree, EdgeInfo **EIs,
  unsigned *ptruNode1, unsigned *ptruNode2,
  double *ptrdLength1, double *ptrdLength2)
	{
	const unsigned uNodeCount = tree.GetNodeCount();

	unsigned uLeaf1 = NULL_NEIGHBOR;
	unsigned uMostDistantLeaf = NULL_NEIGHBOR;
	double dMaxDist = -VERY_LARGE_DOUBLE;
	for (unsigned uNodeIndex = 0; uNodeIndex < uNodeCount; ++uNodeIndex)
		{
		if (!tree.IsLeaf(uNodeIndex))
			continue;

		const unsigned uNode2 = tree.GetNeighbor1(uNodeIndex);
		if (NULL_NEIGHBOR == uNode2)
			Quit("RootByMidLongestSpan: internal error 0");
		const double dEdgeLength = tree.GetEdgeLength(uNodeIndex, uNode2);
		const EdgeInfo &EI = EIs[uNodeIndex][0];
		if (!EI.m_bSet)
			Quit("RootByMidLongestSpan: internal error 1");
		if (EI.m_uNode1 != uNodeIndex || EI.m_uNode2 != uNode2)
			Quit("RootByMidLongestSpan: internal error 2");
		const double dSpanLength = dEdgeLength + EI.m_dMaxDistToLeaf;
		if (dSpanLength > dMaxDist)
			{
			dMaxDist = dSpanLength;
			uLeaf1 = uNodeIndex;
			uMostDistantLeaf = EI.m_uMostDistantLeaf;
			}
		}
	
	if (NULL_NEIGHBOR == uLeaf1)
		Quit("RootByMidLongestSpan: internal error 3");

	const double dTreeHeight = dMaxDist/2.0;
	unsigned uNode1 = uLeaf1;
	unsigned uNode2 = tree.GetNeighbor1(uLeaf1);
	double dAccumSpanLength = 0;

#if	TRACE
	Log("RootByMidLongestSpan: span=%u", uLeaf1);
#endif

	for (;;)
		{
		const double dEdgeLength = tree.GetEdgeLength(uNode1, uNode2);
#if	TRACE
		Log("->%u(%g;%g)", uNode2, dEdgeLength, dAccumSpanLength);
#endif
		if (dAccumSpanLength + dEdgeLength >= dTreeHeight)
			{
			*ptruNode1 = uNode1;
			*ptruNode2 = uNode2;
			*ptrdLength1 = dTreeHeight - dAccumSpanLength;
			*ptrdLength2 = dEdgeLength - *ptrdLength1;
#if	TRACE
			{
			const EdgeInfo &EI = EIs[uLeaf1][0];
			Log("...\n");
			Log("Midpoint: Leaf1=%u Leaf2=%u Node1=%u Node2=%u Length1=%g Length2=%g\n",
			  uLeaf1, EI.m_uMostDistantLeaf, *ptruNode1, *ptruNode2, *ptrdLength1, *ptrdLength2);
			}
#endif
			return;
			}

		if (tree.IsLeaf(uNode2))
			Quit("RootByMidLongestSpan: internal error 4");

		dAccumSpanLength += dEdgeLength;
		const unsigned uSub = tree.GetNeighborSubscript(uNode1, uNode2);
		const EdgeInfo &EI = EIs[uNode1][uSub];
		if (!EI.m_bSet)
			Quit("RootByMidLongestSpan: internal error 5");

		uNode1 = uNode2;
		uNode2 = EI.m_uMaxStep;
		}
	}

/***
Root by balancing average distance to leaves.
The root is a point p such that the average
distance to leaves to the left of p is the
same as the to the right.

This is the method used by CLUSTALW, which
was originally used in PROFILEWEIGHT:

	Thompson et al. (1994) CABIOS (10) 1, 19-29.
***/

static void RootByMinAvgLeafDist(const Tree &tree, EdgeInfo **EIs,
  unsigned *ptruNode1, unsigned *ptruNode2,
  double *ptrdLength1, double *ptrdLength2)
	{
	const unsigned uNodeCount = tree.GetNodeCount();
	const unsigned uLeafCount = tree.GetLeafCount();
	unsigned uNode1 = NULL_NEIGHBOR;
	unsigned uNode2 = NULL_NEIGHBOR;
	double dMinHeight = VERY_LARGE_DOUBLE;
	double dBestLength1 = VERY_LARGE_DOUBLE;
	double dBestLength2 = VERY_LARGE_DOUBLE;

	for (unsigned uNodeIndex = 0; uNodeIndex < uNodeCount; ++uNodeIndex)
		{
		const unsigned uNeighborCount = tree.GetNeighborCount(uNodeIndex);
		for (unsigned uSub = 0; uSub < uNeighborCount; ++uSub)
			{
			const unsigned uNeighborIndex = tree.GetNeighbor(uNodeIndex, uSub);

		// Avoid visiting same edge a second time in reversed order.
			if (uNeighborIndex < uNodeIndex)
				continue;

			const unsigned uSubRev = tree.GetNeighborSubscript(uNeighborIndex, uNodeIndex);
			if (NULL_NEIGHBOR == uSubRev)
				Quit("RootByMinAvgLeafDist, internal error 1");

		// Get info for edges Node1->Node2 and Node2->Node1 (reversed)
			const EdgeInfo &EI = EIs[uNodeIndex][uSub];
			const EdgeInfo &EIRev = EIs[uNeighborIndex][uSubRev];

			if (EI.m_uNode1 != uNodeIndex || EI.m_uNode2 != uNeighborIndex ||
			  EIRev.m_uNode1 != uNeighborIndex || EIRev.m_uNode2 != uNodeIndex)
				Quit("RootByMinAvgLeafDist, internal error 2");
			if (!EI.m_bSet)
				Quit("RootByMinAvgLeafDist, internal error 3");
			if (uLeafCount != EI.m_uLeafCount + EIRev.m_uLeafCount)
				Quit("RootByMinAvgLeafDist, internal error 4");

			const double dEdgeLength = tree.GetEdgeLength(uNodeIndex, uNeighborIndex);
			if (dEdgeLength != tree.GetEdgeLength(uNeighborIndex, uNodeIndex))
				Quit("RootByMinAvgLeafDist, internal error 5");

		// Consider point p on edge 12 in tree (1=Node, 2=Neighbor).
		//
        //	-----         ----
        //	     |       |
        //	     1----p--2
        //	     |       |
        //	-----         ----
		//
		// Define:
		//    ADLp = average distance to leaves to left of point p.
		//	  ADRp = average distance to leaves to right of point p.
		//	  L = edge length = distance 12
		//    x = distance 1p
		// So distance p2 = L - x.
		// Average distance from p to leaves on left of p is:
		//		ADLp = ADL1 + x
		// Average distance from p to leaves on right of p is:
		//		ADRp = ADR2 + (L - x)
		// To be a root, we require these two distances to be equal,
		//		ADLp = ADRp
		//		ADL1 + x = ADR2 + (L - x)
		// Solving for x,
		//		x = (ADR2 - ADL1 + L)/2
		// If 0 <= x <= L, we can place the root on edge 12.

			const double ADL1 = EI.m_dTotalDistToLeaves / EI.m_uLeafCount;
			const double ADR2 = EIRev.m_dTotalDistToLeaves / EIRev.m_uLeafCount;

			const double x = (ADR2 - ADL1 + dEdgeLength)/2.0;
			if (x >= 0 && x <= dEdgeLength)
				{
				const double dLength1 = x;
				const double dLength2 = dEdgeLength - x;
				const double dHeight1 = EI.m_dMaxDistToLeaf + dLength1;
				const double dHeight2 = EIRev.m_dMaxDistToLeaf + dLength2;
				const double dHeight = dHeight1 >= dHeight2 ? dHeight1 : dHeight2;
#if	TRACE
				Log("Candidate root Node1=%u Node2=%u Height=%g\n",
				  uNodeIndex, uNeighborIndex, dHeight);
#endif
				if (dHeight < dMinHeight)
					{
					uNode1 = uNodeIndex;
					uNode2 = uNeighborIndex;
					dBestLength1 = dLength1;
					dBestLength2 = dLength2;
					dMinHeight = dHeight;
					}
				}
			}
		}

	if (NULL_NEIGHBOR == uNode1 || NULL_NEIGHBOR == uNode2)
		Quit("RootByMinAvgLeafDist, internal error 6");

#if	TRACE
	Log("Best root Node1=%u Node2=%u Length1=%g Length2=%g Height=%g\n",
	  uNode1, uNode2, dBestLength1, dBestLength2, dMinHeight);
#endif

	*ptruNode1 = uNode1;
	*ptruNode2 = uNode2;
	*ptrdLength1 = dBestLength1;
	*ptrdLength2 = dBestLength2;
	}

void FixRoot(Tree &tree, ROOT Method)
	{
	if (!tree.IsRooted())
		Quit("FixRoot: expecting rooted tree");

	// Pseudo-root: keep root assigned by clustering
	if (ROOT_FromClustering == Method)
		return;

	tree.UnrootByDeletingRoot();
	tree.RootUnrootedTree(Method);
	}
