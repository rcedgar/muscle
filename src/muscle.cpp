#include "muscle.h"
#include "msa.h"
#include "seqvect.h"
#include "msa.h"
#include "tree.h"
#include "profile.h"

void MUSCLE(SeqVect &v, MSA &msaOut)
	{
	const unsigned uSeqCount = v.Length();

	if (0 == uSeqCount)
		Quit("No sequences in input file");

	ALPHA Alpha = ALPHA_Undefined;
	switch (g_SeqType)
		{
	case SEQTYPE_Auto:
		Alpha = v.GuessAlpha();
		break;

	case SEQTYPE_Protein:
		Alpha = ALPHA_Amino;
		break;

	case SEQTYPE_RNA:
		Alpha = ALPHA_RNA;
		break;

	case SEQTYPE_DNA:
		Alpha = ALPHA_DNA;
		break;

	default:
		Quit("Invalid seq type");
		}
	SetAlpha(Alpha);
	v.FixAlpha();

	if (ALPHA_DNA == Alpha || ALPHA_RNA == Alpha)
		{
		SetPPScore(PPSCORE_SPN);
		g_Distance1 = DISTANCE_Kmer4_6;
		}

	unsigned uMinL = 0;
	unsigned uMaxL = 0;
	unsigned uTotL = 0;
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		unsigned L = v.GetSeq(uSeqIndex).Length();
		uTotL += L;
		if (uMinL == 0 || L < uMinL)
			uMinL = L;
		if (L > uMaxL)
			uMaxL = L;
		}

	SetIter(1);
	g_bDiags = g_bDiags1;
	ShowSeqStats(uSeqCount, uMinL, uMaxL, uTotL/uSeqCount);

	MSA::SetIdCount(uSeqCount);

//// Initialize sequence ids.
//// From this point on, ids must somehow propogate from here.
//	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
//		v.SetSeqId(uSeqIndex, uSeqIndex);

	if (uSeqCount > 1)
		MHackStart(v);

	if (0 == uSeqCount)
		{
		msaOut.Clear();
		return;
		}

	if (1 == uSeqCount && ALPHA_Amino == Alpha)
		{
		const Seq &s = v.GetSeq(0);
		msaOut.FromSeq(s);
		return;
		}

// First iteration
	Tree GuideTree;
	TreeFromSeqVect(v, GuideTree, g_Cluster1, g_Distance1, g_Root1);

	SetMuscleTree(GuideTree);

	ProgNode *ProgNodes = 0;
	ProgressiveAlign(v, GuideTree, msaOut);
	SetCurrentAlignment(msaOut);

	if (1 == g_uMaxIters || 2 == uSeqCount)
		{
		MHackEnd(msaOut);
		return;
		}

	g_bDiags = g_bDiags2;
	SetIter(2);

	RefineTree(msaOut, GuideTree);

	extern void DeleteProgNode(ProgNode &Node);
	const unsigned uNodeCount = GuideTree.GetNodeCount();
	for (unsigned uNodeIndex = 0; uNodeIndex < uNodeCount; ++uNodeIndex)
		DeleteProgNode(ProgNodes[uNodeIndex]);

	delete[] ProgNodes;
	ProgNodes = 0;

	SetSeqWeightMethod(g_SeqWeight2);
	SetMuscleTree(GuideTree);

	if (g_bAnchors)
		RefineVert(msaOut, GuideTree, g_uMaxIters - 2);
	else
		RefineHoriz(msaOut, GuideTree, g_uMaxIters - 2, false, false);

	MHackEnd(msaOut);
	}
