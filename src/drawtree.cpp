#include "muscle.h"
#include "tree.h"

/***
Simple tree drawing algorithm.

y coordinate of node is index in depth-first traversal.
x coordinate is distance from root.
***/

static unsigned DistFromRoot(const Tree &tree, unsigned uNodeIndex)
	{
	const unsigned uRoot = tree.GetRootNodeIndex();
	unsigned uDist = 0;
	while (uNodeIndex != uRoot)
		{
		++uDist;
		uNodeIndex = tree.GetParent(uNodeIndex);
		}
	return uDist;
	}

static void DrawNode(const Tree &tree, unsigned uNodeIndex)
	{
	if (!tree.IsLeaf(uNodeIndex))
		DrawNode(tree, tree.GetLeft(uNodeIndex));

	unsigned uDist = DistFromRoot(tree, uNodeIndex);
	for (unsigned i = 0; i < 5*uDist; ++i)
		Log(" ");
	Log("%d\n", uNodeIndex);

	if (!tree.IsLeaf(uNodeIndex))
		DrawNode(tree, tree.GetRight(uNodeIndex));
	}

void DrawTree(const Tree &tree)
	{
	unsigned uRoot = tree.GetRootNodeIndex();
	DrawNode(tree, uRoot);
	}
