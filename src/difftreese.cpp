#include "muscle.h"
#include "tree.h"

#define TRACE	0

/***
Algorithm to compare two trees, X and Y.

A node x in X and node y in Y are defined to be
similar iff the set of leaves in the subtree under
x is identical to the set of leaves under y.

A node is defined to be changed iff it is not
similar to any node in the other tree.

Nodes x and y are defined to be married iff every
node in the subtree under x is similar to a node
in the subtree under y. Married nodes are considered
to be equal. The subtrees under two married nodes can
at most differ by exchanges of left and right branches,
which we do not consider to be significant here.

A node is changed iff it is not married. If a node is
changed, then it has a dissimilar node in its subtree,
and it follows immediately from the definition of marriage
that its parent is also a bachelor. Hence all nodes on the
path from a changed node to the root are changed.

We assume the trees have the same set of leaves, so
every leaf is trivially both similar and married to
the same leaf in the opposite tree. Changed nodes
are therefore always internal (i.e., non-leaf) nodes.

Example:

              -----A
        -----k
   ----j      -----B
--i     -----C
   ------D


              -----A
        -----p
   ----n      -----B
--m     -----D
   ------C


The following pairs of internal nodes are similar.

	Nodes	Set of leaves
	-----	-------------
	k,p		A,B
	i,m		A,B,C,D

Changed nodes in the first tree are i and j, changed nodes
in the second tree are m and n.

Node k and p are married, but i and m are not (because j
and n are changed). The diffs are C, D and k.

To achieve O(N) we avoid traversing a given subtree multiple
times and also avoid comparing lists of leaves. 

We visit nodes in depth-first order (i.e., a node is visited
before its parent).

If either child of a node is changed, we flag it as changed.

If both children of the node we are visiting are married,
we check whether the spouses of those children have the
same parent in the other tree. If the parents are different,
the current node is a bachelor. If they have the same parent,
then the node we are visiting is the spouse of that parent.
We assign this newly identified married couple a unique integer
id. The id of a node is in one-to-one correspondence with the
set of leaves in its subtree. Two nodes have the same set of
leaves iff they have the same id. Changed nodes do not get
an id.
***/

void DiffTreesE(const Tree &NewTree, const Tree &OldTree,
  unsigned NewNodeIndexToOldNodeIndex[])
	{
#if	TRACE
	Log("DiffTreesE NewTree:\n");
	NewTree.LogMe();
	Log("\n");
	Log("OldTree:\n");
	OldTree.LogMe();
#endif

	if (!NewTree.IsRooted() || !OldTree.IsRooted())
		Quit("DiffTrees: requires rooted trees");

	const unsigned uNodeCount = NewTree.GetNodeCount();
	const unsigned uOldNodeCount = OldTree.GetNodeCount();
	const unsigned uLeafCount = NewTree.GetLeafCount();
	const unsigned uOldLeafCount = OldTree.GetLeafCount();
	if (uNodeCount != uOldNodeCount || uLeafCount != uOldLeafCount)
		Quit("DiffTreesE: different node counts");

	{
	unsigned *IdToOldNodeIndex = new unsigned[uNodeCount];
	for (unsigned uOldNodeIndex = 0; uOldNodeIndex < uNodeCount; ++uOldNodeIndex)
		{
		if (OldTree.IsLeaf(uOldNodeIndex))
			{
			unsigned Id = OldTree.GetLeafId(uOldNodeIndex);
			IdToOldNodeIndex[Id] = uOldNodeIndex;
			}
		}

// Initialize NewNodeIndexToOldNodeIndex[]
// All internal nodes are marked as changed, but may be updated later.
	for (unsigned uNewNodeIndex = 0; uNewNodeIndex < uNodeCount; ++uNewNodeIndex)
		{
		if (NewTree.IsLeaf(uNewNodeIndex))
			{
			unsigned uId = NewTree.GetLeafId(uNewNodeIndex);
			assert(uId < uLeafCount);

			unsigned uOldNodeIndex = IdToOldNodeIndex[uId];
			assert(uOldNodeIndex < uNodeCount);

			NewNodeIndexToOldNodeIndex[uNewNodeIndex] = uOldNodeIndex;
			}
		else
			NewNodeIndexToOldNodeIndex[uNewNodeIndex] = NODE_CHANGED;
		}
	delete[] IdToOldNodeIndex;
	}

// Depth-first traversal of tree.
// The order guarantees that a node is visited before
// its parent is visited.
	for (unsigned uNewNodeIndex = NewTree.FirstDepthFirstNode();
	  NULL_NEIGHBOR != uNewNodeIndex;
	  uNewNodeIndex = NewTree.NextDepthFirstNode(uNewNodeIndex))
		{
		if (NewTree.IsLeaf(uNewNodeIndex))
			continue;

	// If either child is changed, flag this node as changed and continue.
		unsigned uNewLeft = NewTree.GetLeft(uNewNodeIndex);
		unsigned uOldLeft = NewNodeIndexToOldNodeIndex[uNewLeft];
		if (NODE_CHANGED == uOldLeft)
			{
			NewNodeIndexToOldNodeIndex[uNewLeft] = NODE_CHANGED;
			continue;
			}

		unsigned uNewRight = NewTree.GetRight(uNewNodeIndex);
		unsigned uOldRight = NewNodeIndexToOldNodeIndex[uNewRight];
		if (NODE_CHANGED == NewNodeIndexToOldNodeIndex[uNewRight])
			{
			NewNodeIndexToOldNodeIndex[uNewRight] = NODE_CHANGED;
			continue;
			}

		unsigned uOldParentLeft = OldTree.GetParent(uOldLeft);
		unsigned uOldParentRight = OldTree.GetParent(uOldRight);
		if (uOldParentLeft == uOldParentRight)
			NewNodeIndexToOldNodeIndex[uNewNodeIndex] = uOldParentLeft;
		else
			NewNodeIndexToOldNodeIndex[uNewNodeIndex] = NODE_CHANGED;
		}

#if TRACE
	{
	Log("NewToOld ");
	for (unsigned uNewNodeIndex = 0; uNewNodeIndex < uNodeCount; ++uNewNodeIndex)
		{
		Log(" [%3u]=", uNewNodeIndex);
		if (NODE_CHANGED == NewNodeIndexToOldNodeIndex[uNewNodeIndex])
			Log("  X");
		else
			Log("%3u", NewNodeIndexToOldNodeIndex[uNewNodeIndex]);
		if ((uNewNodeIndex+1)%8 == 0)
			Log("\n         ");
		}
	Log("\n");
	}
#endif

#if	DEBUG
	{
	for (unsigned uNewNodeIndex = 0; uNewNodeIndex < uNodeCount; ++uNewNodeIndex)
		{
		unsigned uOld = NewNodeIndexToOldNodeIndex[uNewNodeIndex];
		if (NewTree.IsLeaf(uNewNodeIndex))
			{
			if (uOld >= uNodeCount)
				{
				Log("NewNode=%u uOld=%u > uNodeCount=%u\n",
				  uNewNodeIndex, uOld, uNodeCount);
				Quit("Diff check failed");
				}
			unsigned uIdNew = NewTree.GetLeafId(uNewNodeIndex);
			unsigned uIdOld = OldTree.GetLeafId(uOld);
			if (uIdNew != uIdOld)
				{
				Log("NewNode=%u uOld=%u IdNew=%u IdOld=%u\n",
				  uNewNodeIndex, uOld, uIdNew, uIdOld);
				Quit("Diff check failed");
				}
			continue;
			}

		if (NODE_CHANGED == uOld)
			continue;

		unsigned uNewLeft = NewTree.GetLeft(uNewNodeIndex);
		unsigned uNewRight = NewTree.GetRight(uNewNodeIndex);

		unsigned uOldLeft = OldTree.GetLeft(uOld);
		unsigned uOldRight = OldTree.GetRight(uOld);

		unsigned uNewLeftPartner = NewNodeIndexToOldNodeIndex[uNewLeft];
		unsigned uNewRightPartner = NewNodeIndexToOldNodeIndex[uNewRight];

		bool bSameNotRotated = (uNewLeftPartner == uOldLeft && uNewRightPartner == uOldRight);
		bool bSameRotated = (uNewLeftPartner == uOldRight && uNewRightPartner == uOldLeft);
		if (!bSameNotRotated && !bSameRotated)
			{
			Log("NewNode=%u NewL=%u NewR=%u\n", uNewNodeIndex, uNewLeft, uNewRight);
			Log("OldNode=%u OldL=%u OldR=%u\n", uOld, uOldLeft, uOldRight);
			Log("NewLPartner=%u NewRPartner=%u\n", uNewLeftPartner, uNewRightPartner);
			Quit("Diff check failed");
			}
		}
	}
#endif
	}
