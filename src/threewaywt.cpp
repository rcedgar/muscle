#include "muscle.h"
#include "tree.h"
#include <math.h>

#define TRACE	0

/***
Sequence weights derived from a tree using Gotoh's
three-way method.

	Gotoh (1995) CABIOS 11(5), 543-51.

Each edge e is assigned a weight w(e).

Consider first a tree with three leaves A,B and C
having branch lengths a, b and c, as follows.

            B
            |
            b
            |
    A---a---R---c---C

The internal node is denoted by R.

Define:

	S = (ab + ca + ab)
	x = bc(a + b)(a + c)
	y = a(b + c)FS

Here F is a tunable normalization factor which is
approximately 1.0. Then the edge weight for AR
is computed as:

	w(AR) = sqrt(x/y)

Similar expressions for the other edges follow by
symmetry.

For a tree with more than three edges, the weight
of an edge that ends in a leaf is computed from
the three-way tree that includes the edge and
its two neighbors. The weight of an internal edge
is computed as the product of the weights for that
edge derived from the two three-way subtrees that
include that edge.

For example, consider the following tree.

       B
       |
    A--R--V--C
          |
          D

Here, w(RV) is computed as the product of the
two values for w(RV) derived from the three-way
trees with leaves ABV and RCD respectively.

The calculation is done using "Gotoh lengths",
not the real edge lengths.

The Gotoh length G of a directed edge is calculated
recursively as:

	G = d + LR/(L + R)

where d is the length of the edge, and L and R are
the Gotoh lengths of the left and right edges adjoining
the terminal end of the edge. If the edge terminates on
a leaf, then G=d.

Pairwise sequence weights are computed as the
product of edge weights on the path that connects
their leaves.

If the tree is split into two subtrees by deleting
a given edge e, then the pairwise weights factorize.
For operations on profiles formed from the two
subtrees, it is possible to assign a weight to a
sequence as the product of edge weights on a path
from e to its leaf.
***/

// The xxxUnrooted functions present a rooted tree as
// if it had been unrooted by deleting the root node.
static unsigned GetFirstNeighborUnrooted(const Tree &tree, unsigned uNode1,
  unsigned uNode2)
	{
	if (tree.IsRoot(uNode1) || tree.IsRoot(uNode2))
		Quit("GetFirstNeighborUnrooted, should never be called with root");
	if (!tree.IsEdge(uNode1, uNode2))
		{
		if (!tree.IsRoot(tree.GetParent(uNode1)) ||
		  !tree.IsRoot(tree.GetParent(uNode2)))
			Quit("GetFirstNeighborUnrooted, not edge");
		const unsigned uRoot = tree.GetRootNodeIndex();
		return tree.GetFirstNeighbor(uNode1, uRoot);
		}

	unsigned uNeighbor = tree.GetFirstNeighbor(uNode1, uNode2);
	if (tree.IsRoot(uNeighbor))
		return tree.GetFirstNeighbor(uNeighbor, uNode1);
	return uNeighbor;
	}

static unsigned GetSecondNeighborUnrooted(const Tree &tree, unsigned uNode1,
  unsigned uNode2)
	{
	if (tree.IsRoot(uNode1) || tree.IsRoot(uNode2))
		Quit("GetFirstNeighborUnrooted, should never be called with root");
	if (!tree.IsEdge(uNode1, uNode2))
		{
		if (!tree.IsRoot(tree.GetParent(uNode1)) ||
		  !tree.IsRoot(tree.GetParent(uNode2)))
			Quit("GetFirstNeighborUnrooted, not edge");
		const unsigned uRoot = tree.GetRootNodeIndex();
		return tree.GetSecondNeighbor(uNode1, uRoot);
		}

	unsigned uNeighbor = tree.GetSecondNeighbor(uNode1, uNode2);
	if (tree.IsRoot(uNeighbor))
		return tree.GetFirstNeighbor(uNeighbor, uNode1);
	return uNeighbor;
	}

static unsigned GetNeighborUnrooted(const Tree &tree, unsigned uNode1,
  unsigned uSub)
	{
	unsigned uNeighbor = tree.GetNeighbor(uNode1, uSub);
	if (tree.IsRoot(uNeighbor))
		return tree.GetFirstNeighbor(uNeighbor, uNode1);
	return uNeighbor;
	}

static unsigned GetNeighborSubscriptUnrooted(const Tree &tree, unsigned uNode1,
  unsigned uNode2)
	{
	if (tree.IsEdge(uNode1, uNode2))
		return tree.GetNeighborSubscript(uNode1, uNode2);
	if (!tree.IsRoot(tree.GetParent(uNode1)) ||
	  !tree.IsRoot(tree.GetParent(uNode2)))
		Quit("GetNeighborSubscriptUnrooted, not edge");
	for (unsigned uSub = 0; uSub < 3; ++uSub)
		if (GetNeighborUnrooted(tree, uNode1, uSub) == uNode2)
			return uSub;
	Quit("GetNeighborSubscriptUnrooted, not a neighbor");
	return NULL_NEIGHBOR;
	}

static double GetEdgeLengthUnrooted(const Tree &tree, unsigned uNode1,
  unsigned uNode2)
	{
	if (tree.IsRoot(uNode1) || tree.IsRoot(uNode2))
		Quit("GetEdgeLengthUnrooted, should never be called with root");
	if (!tree.IsEdge(uNode1, uNode2))
		{
		if (!tree.IsRoot(tree.GetParent(uNode1)) ||
		  !tree.IsRoot(tree.GetParent(uNode2)))
			Quit("GetEdgeLengthUnrooted, not edge");

		const unsigned uRoot = tree.GetRootNodeIndex();
		return tree.GetEdgeLength(uNode1, uRoot) +
		  tree.GetEdgeLength(uNode2, uRoot);
		}
	return tree.GetEdgeLength(uNode1, uNode2);
	}

double GetGotohLength(const Tree &tree, unsigned R, unsigned A)
	{
	double dThis = GetEdgeLengthUnrooted(tree, R, A);

// Enforce non-negative edge lengths
	if (dThis < 0)
		dThis = 0;

	if (tree.IsLeaf(A))
		return dThis;

	const unsigned uFirst = GetFirstNeighborUnrooted(tree, A, R);
	const unsigned uSecond = GetSecondNeighborUnrooted(tree, A, R);
	const double dFirst = GetGotohLength(tree, A, uFirst);
	const double dSecond = GetGotohLength(tree, A, uSecond);
	const double dSum = dFirst + dSecond;
	const double dThird = dSum == 0 ? 0 : (dFirst*dSecond)/dSum;
	return dThis + dThird;
	}

// Return weight of edge A-R in three-way subtree that has
// leaves A,B,C and internal node R.
static double GotohWeightThreeWay(const Tree &tree, unsigned A,
  unsigned B, unsigned C, unsigned R)
	{
	const double F = 1.0;

	if (tree.IsLeaf(R))
		Quit("GotohThreeWay: R must be internal node");

	double a = GetGotohLength(tree, R, A);
	double b = GetGotohLength(tree, R, B);
	double c = GetGotohLength(tree, R, C);

	double S = b*c + c*a + a*b;
	double x = b*c*(a + b)*(a + c);
	double y = a*(b + c)*F*S;

// y is zero iff all three branch lengths are zero.
	if (y < 0.001)
		return 1.0;
	return sqrt(x/y);
	}

static double GotohWeightEdge(const Tree &tree, unsigned uNodeIndex1,
  unsigned uNodeIndex2)
	{
	double w1 = 1.0;
	double w2 = 1.0;
	if (!tree.IsLeaf(uNodeIndex1))
		{
		unsigned R = uNodeIndex1;
		unsigned A = uNodeIndex2;
		unsigned B = GetFirstNeighborUnrooted(tree, R, A);
		unsigned C = GetSecondNeighborUnrooted(tree, R, A);
		w1 = GotohWeightThreeWay(tree, A, B, C, R);
		}
	if (!tree.IsLeaf(uNodeIndex2))
		{
		unsigned R = uNodeIndex2;
		unsigned A = uNodeIndex1;
		unsigned B = GetFirstNeighborUnrooted(tree, R, A);
		unsigned C = GetSecondNeighborUnrooted(tree, R, A);
		w2 = GotohWeightThreeWay(tree, A, B, C, R);
		}
	return w1*w2;
	}

void CalcThreeWayEdgeWeights(const Tree &tree, WEIGHT **EdgeWeights)
	{
	const unsigned uNodeCount = tree.GetNodeCount();
	for (unsigned uNodeIndex1 = 0; uNodeIndex1 < uNodeCount; ++uNodeIndex1)
		{
		if (tree.IsRoot(uNodeIndex1))
			continue;
		for (unsigned uSub1 = 0; uSub1 < 3; ++uSub1)
			{
			const unsigned uNodeIndex2 = GetNeighborUnrooted(tree, uNodeIndex1, uSub1);
			if (NULL_NEIGHBOR == uNodeIndex2)
				continue;

		// Avoid computing same edge twice in reversed order
			if (uNodeIndex2 < uNodeIndex1)
				continue;

			const WEIGHT w = (WEIGHT) GotohWeightEdge(tree, uNodeIndex1, uNodeIndex2);
			const unsigned uSub2 = GetNeighborSubscriptUnrooted(tree, uNodeIndex2, uNodeIndex1);
#if	DEBUG
			{
			assert(uNodeIndex2 == GetNeighborUnrooted(tree, uNodeIndex1, uSub1));
			assert(uNodeIndex1 == GetNeighborUnrooted(tree, uNodeIndex2, uSub2));
			const WEIGHT wRev = (WEIGHT) GotohWeightEdge(tree, uNodeIndex2, uNodeIndex1);
			if (!BTEq(w, wRev))
				Quit("CalcThreeWayWeights: rev check failed %g %g",
				  w, wRev);
			}
#endif
			EdgeWeights[uNodeIndex1][uSub1] = w;
			EdgeWeights[uNodeIndex2][uSub2] = w;
			}
		}
	}

static void SetSeqWeights(const Tree &tree, unsigned uNode1, unsigned uNode2,
  double dPathWeight, WEIGHT *Weights)
	{
	if (tree.IsRoot(uNode1) || tree.IsRoot(uNode2))
		Quit("SetSeqWeights, should never be called with root");

	const double dThisLength = GetEdgeLengthUnrooted(tree, uNode1, uNode2);
	if (tree.IsLeaf(uNode2))
		{
		const unsigned Id = tree.GetLeafId(uNode2);
		Weights[Id] = (WEIGHT) (dPathWeight + dThisLength);
		return;
		}
	const unsigned uFirst = GetFirstNeighborUnrooted(tree, uNode2, uNode1);
	const unsigned uSecond = GetSecondNeighborUnrooted(tree, uNode2, uNode1);
	dPathWeight *= dThisLength;
	SetSeqWeights(tree, uNode2, uFirst, dPathWeight, Weights);
	SetSeqWeights(tree, uNode2, uSecond, dPathWeight, Weights);
	}

void CalcThreeWayWeights(const Tree &tree, unsigned uNode1, unsigned uNode2,
  WEIGHT *Weights)
	{
#if	TRACE
	Log("CalcThreeWayEdgeWeights\n");
	tree.LogMe();
#endif

	if (tree.IsRoot(uNode1))
		uNode1 = tree.GetFirstNeighbor(uNode1, uNode2);
	else if (tree.IsRoot(uNode2))
		uNode2 = tree.GetFirstNeighbor(uNode2, uNode1);
	const unsigned uNodeCount = tree.GetNodeCount();
	WEIGHT **EdgeWeights = new WEIGHT *[uNodeCount];
	for (unsigned uNodeIndex = 0; uNodeIndex < uNodeCount; ++uNodeIndex)
		EdgeWeights[uNodeIndex] = new WEIGHT[3];

	CalcThreeWayEdgeWeights(tree, EdgeWeights);

#if	TRACE
	{
	Log("Node1  Node2  Length   Gotoh  EdgeWt\n");
	Log("-----  -----  ------  ------  ------\n");
	for (unsigned uNodeIndex1 = 0; uNodeIndex1 < uNodeCount; ++uNodeIndex1)
		{
		if (tree.IsRoot(uNodeIndex1))
			continue;
		for (unsigned uSub1 = 0; uSub1 < 3; ++uSub1)
			{
			const unsigned uNodeIndex2 = GetNeighborUnrooted(tree, uNodeIndex1, uSub1);
			if (NULL_NEIGHBOR == uNodeIndex2)
				continue;
			if (uNodeIndex2 < uNodeIndex1)
				continue;
			const WEIGHT ew = EdgeWeights[uNodeIndex1][uSub1];
			const double d = GetEdgeLengthUnrooted(tree, uNodeIndex1, uNodeIndex2);
			const double g = GetGotohLength(tree, uNodeIndex1, uNodeIndex2);
			Log("%5u  %5u  %6.3f  %6.3f  %6.3f\n", uNodeIndex1, uNodeIndex2, d, g, ew);
			}
		}
	}
#endif

	SetSeqWeights(tree, uNode1, uNode2, 0.0, Weights);
	SetSeqWeights(tree, uNode2, uNode1, 0.0, Weights);

	for (unsigned uNodeIndex = 0; uNodeIndex < uNodeCount; ++uNodeIndex)
		delete[] EdgeWeights[uNodeIndex];
	delete[] EdgeWeights;
	}
