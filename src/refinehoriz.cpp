#include "muscle.h"
#include "tree.h"
#include "msa.h"
#include "pwpath.h"
#include "profile.h"
#include "scorehistory.h"
#include "objscore.h"

unsigned g_uRefineHeightSubtree;
unsigned g_uRefineHeightSubtreeTotal;

#define TRACE			0
#define DIFFOBJSCORE	0

static bool TryRealign(MSA &msaIn, const Tree &tree, const unsigned Leaves1[],
  unsigned uCount1, const unsigned Leaves2[], unsigned uCount2,
  SCORE *ptrscoreBefore, SCORE *ptrscoreAfter,
  bool bLockLeft, bool bLockRight)
	{
#if	TRACE
	Log("TryRealign, msaIn=\n");
	msaIn.LogMe();
#endif

	const unsigned uSeqCount = msaIn.GetSeqCount();

	unsigned *Ids1 = new unsigned[uSeqCount];
	unsigned *Ids2 = new unsigned[uSeqCount];

	LeafIndexesToIds(tree, Leaves1, uCount1, Ids1);
	LeafIndexesToIds(tree, Leaves2, uCount2, Ids2);

	MSA msa1;
	MSA msa2;

	MSASubsetByIds(msaIn, Ids1, uCount1, msa1);
	MSASubsetByIds(msaIn, Ids2, uCount2, msa2);

#if	DEBUG
	ValidateMuscleIds(msa1);
	ValidateMuscleIds(msa2);
#endif

// Computing the objective score may be expensive for
// large numbers of sequences. As a speed optimization,
// we check whether the alignment changes. If it does
// not change, there is no need to compute the objective
// score. We test for the alignment changing by comparing
// the Viterbi paths before and after re-aligning.
	PWPath pathBefore;
	pathBefore.FromMSAPair(msa1, msa2);

	DeleteGappedCols(msa1);
	DeleteGappedCols(msa2);

	if (0 == msa1.GetColCount() || 0 == msa2.GetColCount())
		return false;

	MSA msaRealigned;
	PWPath pathAfter;

	AlignTwoMSAs(msa1, msa2, msaRealigned, pathAfter, bLockLeft, bLockRight);

	bool bAnyChanges = !pathAfter.Equal(pathBefore);
	unsigned uDiffCount1;
	unsigned uDiffCount2;
	uint EdgeCount1 = pathBefore.GetEdgeCount();
	uint EdgeCount2 = pathAfter.GetEdgeCount();
	uint *Edges1 = new uint[EdgeCount1];
	uint *Edges2 = new uint[EdgeCount2];
	DiffPaths(pathBefore, pathAfter, Edges1, &uDiffCount1, Edges2, &uDiffCount2);
	delete[] Edges1;
	delete[] Edges2;
	Edges1 = 0;
	Edges2 = 0;

#if	TRACE
	Log("TryRealign, msa1=\n");
	msa1.LogMe();
	Log("\nmsa2=\n");
	msa2.LogMe();
	Log("\nRealigned (changes %s)=\n", bAnyChanges ? "TRUE" : "FALSE");
	msaRealigned.LogMe();
#endif

	if (!bAnyChanges)
		{
		*ptrscoreBefore = 0;
		*ptrscoreAfter = 0;
		return false;
		}

	SetMSAWeightsMuscle(msaIn);
	SetMSAWeightsMuscle(msaRealigned);

#if	DIFFOBJSCORE
	const SCORE scoreDiff = DiffObjScore(msaIn, pathBefore, Edges1, uDiffCount1,
	  msaRealigned, pathAfter, Edges2, uDiffCount2);
	bool bAccept = (scoreDiff > 0);
	*ptrscoreBefore = 0;
	*ptrscoreAfter = scoreDiff;
	//const SCORE scoreBefore = ObjScoreIds(msaIn, Ids1, uCount1, Ids2, uCount2);
	//const SCORE scoreAfter = ObjScoreIds(msaRealigned, Ids1, uCount1, Ids2, uCount2);
	//Log("Diff = %.3g %.3g\n", scoreDiff, scoreAfter - scoreBefore);
#else
	const SCORE scoreBefore = ObjScoreIds(msaIn, Ids1, uCount1, Ids2, uCount2);
	const SCORE scoreAfter = ObjScoreIds(msaRealigned, Ids1, uCount1, Ids2, uCount2);

	bool bAccept = (scoreAfter > scoreBefore);

#if	TRACE
	Log("Score %g -> %g Accept %s\n", scoreBefore, scoreAfter, bAccept ? "TRUE" : "FALSE");
#endif

	*ptrscoreBefore = scoreBefore;
	*ptrscoreAfter = scoreAfter;
#endif

	if (bAccept)
		msaIn.Copy(msaRealigned);
	delete[] Ids1;
	delete[] Ids2;
	return bAccept;
	}

static void RefineHeightParts(MSA &msaIn, const Tree &tree,
 const unsigned InternalNodeIndexes[], bool bReversed, bool bRight,
 unsigned uIter, 
 ScoreHistory &History,
 bool *ptrbAnyChanges, bool *ptrbOscillating, bool bLockLeft, bool bLockRight)
	{
	*ptrbOscillating = false;

	const unsigned uSeqCount = msaIn.GetSeqCount();
	const unsigned uInternalNodeCount = uSeqCount - 1;

	unsigned *Leaves1 = new unsigned[uSeqCount];
	unsigned *Leaves2 = new unsigned[uSeqCount];

	const unsigned uRootNodeIndex = tree.GetRootNodeIndex();
	bool bAnyAccepted = false;
	for (unsigned i = 0; i < uInternalNodeCount; ++i)
		{
		const unsigned uInternalNodeIndex = InternalNodeIndexes[i];
		unsigned uNeighborNodeIndex;
		if (tree.IsRoot(uInternalNodeIndex) && !bRight)
			continue;
		else if (bRight)
			uNeighborNodeIndex = tree.GetRight(uInternalNodeIndex);
		else
			uNeighborNodeIndex = tree.GetLeft(uInternalNodeIndex);

		g_uTreeSplitNode1 = uInternalNodeIndex;
		g_uTreeSplitNode2 = uNeighborNodeIndex;

		unsigned uCount1;
		unsigned uCount2;

		GetLeaves(tree, uNeighborNodeIndex, Leaves1, &uCount1);
		GetLeavesExcluding(tree, uRootNodeIndex, uNeighborNodeIndex,
		  Leaves2, &uCount2);

#if	TRACE
		Log("\nRefineHeightParts node %u\n", uInternalNodeIndex);
		Log("Group1=");
		for (unsigned n = 0; n < uCount1; ++n)
			Log(" %u(%s)", Leaves1[n], tree.GetName(Leaves1[n]));
		Log("\n");
		Log("Group2=");
		for (unsigned n = 0; n < uCount2; ++n)
			Log(" %u(%s)", Leaves2[n], tree.GetName(Leaves2[n]));
		Log("\n");
#endif

		SCORE scoreBefore;
		SCORE scoreAfter;
		bool bAccepted = TryRealign(msaIn, tree, Leaves1, uCount1, Leaves2, uCount2,
		  &scoreBefore, &scoreAfter, bLockLeft, bLockRight);
		SetCurrentAlignment(msaIn);

		++g_uRefineHeightSubtree;
		Progress(g_uRefineHeightSubtree, g_uRefineHeightSubtreeTotal);

#if	TRACE
		if (uIter > 0)
			Log("Before %g %g\n", scoreBefore,
			  History.GetScore(uIter - 1, uInternalNodeIndex, bReversed, bRight));
#endif
		SCORE scoreMax = scoreAfter > scoreBefore? scoreAfter : scoreBefore;
		bool bRepeated = History.SetScore(uIter, uInternalNodeIndex, bRight, scoreMax);
		if (bRepeated)
			{
			*ptrbOscillating = true;
			break;
			}

		if (bAccepted)
			bAnyAccepted = true;
		}

	delete[] Leaves1;
	delete[] Leaves2;

	*ptrbAnyChanges = bAnyAccepted;
	}

// Return true if any changes made
bool RefineHoriz(MSA &msaIn, const Tree &tree, unsigned uIters, bool bLockLeft,
  bool bLockRight)
	{
#if	TRACE
	tree.LogMe();
#endif

	if (!tree.IsRooted())
		Quit("RefineHeight: requires rooted tree");

	const unsigned uSeqCount = msaIn.GetSeqCount();
	if (uSeqCount < 3)
		return false;

	const unsigned uInternalNodeCount = uSeqCount - 1;
	unsigned *InternalNodeIndexes = new unsigned[uInternalNodeCount];
	unsigned *InternalNodeIndexesR = new unsigned[uInternalNodeCount];

	GetInternalNodesInHeightOrder(tree, InternalNodeIndexes);

	ScoreHistory History(uIters, 2*uSeqCount - 1);

	bool bAnyChangesAnyIter = false;
	for (unsigned n = 0; n < uInternalNodeCount; ++n)
		InternalNodeIndexesR[uInternalNodeCount - 1 - n] = InternalNodeIndexes[n];

	for (unsigned uIter = 0; uIter < uIters; ++uIter)
		{
		bool bAnyChangesThisIter = false;
		IncIter();
		SetProgressDesc("Refine biparts");
		g_uRefineHeightSubtree = 0;
		g_uRefineHeightSubtreeTotal = uInternalNodeCount*2 - 1;

		bool bReverse = (uIter%2 != 0);
		unsigned *Internals;
		if (bReverse)
			Internals = InternalNodeIndexesR;
		else
			Internals = InternalNodeIndexes;

		bool bOscillating;
		for (unsigned i = 0; i < 2; ++i)
			{
			bool bAnyChanges = false;
			bool bRight;
			switch (i)
				{
			case 0:
				bRight = true;
				break;
			case 1:
				bRight = false;
				break;
			default:
				Quit("RefineHeight default case");
				}
			RefineHeightParts(msaIn, tree, Internals, bReverse, bRight,
			  uIter, 
			  History, 
			  &bAnyChanges, &bOscillating, bLockLeft, bLockRight);
			if (bOscillating)
				{
				ProgressStepsDone();
				goto Osc;
				}
			if (bAnyChanges)
				{
				bAnyChangesThisIter = true;
				bAnyChangesAnyIter = true;
				}
			}

		ProgressStepsDone();
		if (bOscillating)
			break;

		if (!bAnyChangesThisIter)
			break;
		}

Osc:
	delete[] InternalNodeIndexes;
	delete[] InternalNodeIndexesR;

	return bAnyChangesAnyIter;
	}
