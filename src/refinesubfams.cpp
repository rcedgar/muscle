#include "muscle.h"
#include "msa.h"
#include "tree.h"
#include "clust.h"
#include "profile.h"
#include "pwpath.h"

#define TRACE 0

static void ProgressiveAlignSubfams(const Tree &tree, const unsigned Subfams[],
  unsigned uSubfamCount, const MSA SubfamMSAs[], MSA &msa);

// Identify subfamilies in a tree.
// Returns array of internal node indexes, one for each subfamily.
// First try is to select groups by height (which should approximate
// minimum percent identity), if this gives too many subfamilies then
// we cut at a point that gives the maximum allowed number of subfams.
static void GetSubfams(const Tree &tree, double dMaxHeight,
  unsigned uMaxSubfamCount, unsigned **ptrptrSubfams, unsigned *ptruSubfamCount)
	{
	const unsigned uNodeCount = tree.GetNodeCount();

	unsigned *Subfams = new unsigned[uNodeCount];

	unsigned uSubfamCount;
	ClusterByHeight(tree, dMaxHeight, Subfams, &uSubfamCount);

	if (uSubfamCount > uMaxSubfamCount)
		ClusterBySubfamCount(tree, uMaxSubfamCount, Subfams, &uSubfamCount);

	*ptrptrSubfams = Subfams;
	*ptruSubfamCount = uSubfamCount;
	}

static void LogSubfams(const Tree &tree, const unsigned Subfams[],
  unsigned uSubfamCount)
	{
	const unsigned uNodeCount = tree.GetNodeCount();
	Log("%u subfamilies found\n", uSubfamCount);
	Log("Subfam  Sequence\n");
	Log("------  --------\n");
	unsigned *Leaves = new unsigned[uNodeCount];
	for (unsigned uSubfamIndex = 0; uSubfamIndex < uSubfamCount; ++uSubfamIndex)
		{
		unsigned uSubfamNodeIndex = Subfams[uSubfamIndex];
		unsigned uLeafCount;
		GetLeaves(tree, uSubfamNodeIndex, Leaves, &uLeafCount);
		for (unsigned uLeafIndex = 0; uLeafIndex < uLeafCount; ++uLeafIndex)
			Log("%6u  %s\n", uSubfamIndex + 1, tree.GetLeafName(Leaves[uLeafIndex]));
		Log("\n");
		}
	delete[] Leaves;
	}

bool RefineSubfams(MSA &msa, const Tree &tree, unsigned uIters)
	{
	const unsigned uSeqCount = msa.GetSeqCount();
	if (uSeqCount < 3)
		return false;

	const double dMaxHeight = 0.6;
	const unsigned uMaxSubfamCount = 16;
	const unsigned uNodeCount = tree.GetNodeCount();

	unsigned *Subfams;
	unsigned uSubfamCount;
	GetSubfams(tree, dMaxHeight, uMaxSubfamCount, &Subfams, &uSubfamCount);
	assert(uSubfamCount <= uSeqCount);

	if (g_bVerbose)
		LogSubfams(tree, Subfams, uSubfamCount);

	MSA *SubfamMSAs = new MSA[uSubfamCount];
	unsigned *Leaves = new unsigned[uSeqCount];
	unsigned *Ids = new unsigned[uSeqCount];

	bool bAnyChanges = false;
	for (unsigned uSubfamIndex = 0; uSubfamIndex < uSubfamCount; ++uSubfamIndex)
		{
		unsigned uSubfam = Subfams[uSubfamIndex];
		unsigned uLeafCount;
		GetLeaves(tree, uSubfam, Leaves, &uLeafCount);
		assert(uLeafCount <= uSeqCount);

		LeafIndexesToIds(tree, Leaves, uLeafCount, Ids);

		MSA &msaSubfam = SubfamMSAs[uSubfamIndex];
		MSASubsetByIds(msa, Ids, uLeafCount, msaSubfam);
		DeleteGappedCols(msaSubfam);

#if	TRACE
		Log("Subfam %u MSA=\n", uSubfamIndex);
		msaSubfam.LogMe();
#endif

		if (msaSubfam.GetSeqCount() <= 2)
			continue;

	// TODO /////////////////////////////////////////
	// Try using existing tree, may actually hurt to
	// re-estimate, may also be a waste of CPU & mem.
	/////////////////////////////////////////////////
		Tree SubfamTree;
		TreeFromMSA(msaSubfam, SubfamTree, g_Cluster2, g_Distance2, g_Root2);

		bool bAnyChangesThisSubfam;
		if (g_bAnchors)
			bAnyChangesThisSubfam = RefineVert(msaSubfam, SubfamTree, uIters);
		else
			bAnyChangesThisSubfam = RefineHoriz(msaSubfam, SubfamTree, uIters, false, false);
#if	TRACE
		Log("Subfam %u Changed %d\n", uSubfamIndex, bAnyChangesThisSubfam);
#endif
		if (bAnyChangesThisSubfam)
			bAnyChanges = true;
		}

	if (bAnyChanges)
		ProgressiveAlignSubfams(tree, Subfams, uSubfamCount, SubfamMSAs, msa);

	delete[] Leaves;
	delete[] Subfams;
	delete[] SubfamMSAs;

	return bAnyChanges;
	}

static void ProgressiveAlignSubfams(const Tree &tree, const unsigned Subfams[],
  unsigned uSubfamCount, const MSA SubfamMSAs[], MSA &msa)
	{
	const unsigned uNodeCount = tree.GetNodeCount();

	bool *Ready = new bool[uNodeCount];
	MSA **MSAs = new MSA *[uNodeCount];
	for (unsigned uNodeIndex = 0; uNodeIndex < uNodeCount; ++uNodeIndex)
		{
		Ready[uNodeIndex] = false;
		MSAs[uNodeIndex] = 0;
		}

	for (unsigned uSubfamIndex = 0; uSubfamIndex < uSubfamCount; ++uSubfamIndex)
		{
		unsigned uNodeIndex = Subfams[uSubfamIndex];
		Ready[uNodeIndex] = true;
		MSA *ptrMSA = new MSA;
	// TODO: Wasteful copy, needs re-design
		ptrMSA->Copy(SubfamMSAs[uSubfamIndex]);
		MSAs[uNodeIndex] = ptrMSA;
		}

	for (unsigned uNodeIndex = tree.FirstDepthFirstNode();
	  NULL_NEIGHBOR != uNodeIndex;
	  uNodeIndex = tree.NextDepthFirstNode(uNodeIndex))
		{
		if (tree.IsLeaf(uNodeIndex))
			continue;

		unsigned uRight = tree.GetRight(uNodeIndex);
		unsigned uLeft = tree.GetLeft(uNodeIndex);
		if (!Ready[uRight] || !Ready[uLeft])
			continue;

		MSA *ptrLeft = MSAs[uLeft];
		MSA *ptrRight = MSAs[uRight];
		assert(ptrLeft != 0 && ptrRight != 0);

		MSA *ptrParent = new MSA;

		PWPath Path;
		AlignTwoMSAs(*ptrLeft, *ptrRight, *ptrParent, Path);

		MSAs[uNodeIndex] = ptrParent;
		Ready[uNodeIndex] = true;
		Ready[uLeft] = false;
		Ready[uRight] = false;

		delete MSAs[uLeft];
		delete MSAs[uRight];
		MSAs[uLeft] = 0;
		MSAs[uRight] = 0;
		}

#if	DEBUG
	{
	unsigned uReadyCount = 0;
	for (unsigned uNodeIndex = 0; uNodeIndex < uNodeCount; ++uNodeIndex)
		{
		if (Ready[uNodeIndex])
			{
			assert(tree.IsRoot(uNodeIndex));
			++uReadyCount;
			assert(0 != MSAs[uNodeIndex]);
			}
		else
			assert(0 == MSAs[uNodeIndex]);
		}
	assert(1 == uReadyCount);
	}
#endif

	const unsigned uRoot = tree.GetRootNodeIndex();
	MSA *ptrRootAlignment = MSAs[uRoot];

	msa.Copy(*ptrRootAlignment);

	delete ptrRootAlignment;

#if	TRACE
	Log("After refine subfamilies, root alignment=\n");
	msa.LogMe();
#endif
	}
