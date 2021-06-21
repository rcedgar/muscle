#include "muscle.h"
#include "textfile.h"
#include "seqvect.h"
#include "distfunc.h"
#include "msa.h"
#include "tree.h"
#include "profile.h"
#include "timing.h"

static char g_strUseTreeWarning[] =
"\n******** WARNING ****************\n"
"\nYou specified the -usetree option.\n"
"Note that a good evolutionary tree may NOT be a good\n"
"guide tree for multiple alignment. For more details,\n"
"please refer to the user guide. To disable this\n"
"warning, use -usetree_nowarn <treefilename>.\n\n";

void SetAlphaFromSeqVect(SeqVect &v)
	{
	ALPHA Alpha = ALPHA_Undefined;
	switch (g_SeqType)
		{
	case SEQTYPE_Auto:
		Alpha = v.GuessAlpha();
		break;

	case SEQTYPE_Protein:
		Alpha = ALPHA_Amino;
		break;

	case SEQTYPE_DNA:
		Alpha = ALPHA_DNA;
		break;

	case SEQTYPE_RNA:
		Alpha = ALPHA_RNA;
		break;

	default:
		Quit("Invalid seq type");
		}
	SetAlpha(Alpha);
	v.FixAlpha();
	}

void ShowSeqStats(const SeqVect &v)
	{
	const uint SeqCount = v.Length();
	uint uMinL = 0;
	uint uMaxL = 0;
	uint uTotL = 0;
	for (uint uSeqIndex = 0; uSeqIndex < SeqCount; ++uSeqIndex)
		{
		uint L = v.GetSeq(uSeqIndex).Length();
		uTotL += L;
		if (uMinL == 0 || L < uMinL)
			uMinL = L;
		if (L > uMaxL)
			uMaxL = L;
		}
	ShowSeqStats(SeqCount, uMinL, uMaxL, uTotL/SeqCount);
	}

void DoMuscle()
	{
	SetStartTime();
	SetOutputFileName(g_pstrOutFileName);
	SetInputFileName(g_pstrInFileName);

	SetMaxIters(g_uMaxIters);
	SetSeqWeightMethod(g_SeqWeight1);

	TextFile fileIn(g_pstrInFileName);
	SeqVect v;
	v.FromFASTAFile(fileIn);
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

	case SEQTYPE_DNA:
		Alpha = ALPHA_DNA;
		break;

	case SEQTYPE_RNA:
		Alpha = ALPHA_RNA;
		break;

	default:
		Quit("Invalid seq type");
		}
	SetAlpha(Alpha);
	v.FixAlpha();

	PTR_SCOREMATRIX UserMatrix = 0;
	if (0 != g_pstrMatrixFileName)
		{
		const char *FileName = g_pstrMatrixFileName;
		const char *Path = getenv("MUSCLE_MXPATH");
		if (Path != 0)
			{
			size_t n = strlen(Path) + 1 + strlen(FileName) + 1;
			char *NewFileName = new char[n];
			sprintf(NewFileName, "%s/%s", Path, FileName);
			FileName = NewFileName;
			}
		TextFile File(FileName);
		UserMatrix = ReadMx(File);
		g_Alpha = ALPHA_Amino;
		g_PPScore = PPSCORE_SP;
		}

	SetPPScore();

	if (0 != UserMatrix)
		g_ptrScoreMatrix = UserMatrix;

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

	SetMuscleSeqVect(v);

	MSA::SetIdCount(uSeqCount);

// Initialize sequence ids.
// From this point on, ids must somehow propogate from here.
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		v.SetSeqId(uSeqIndex, uSeqIndex);

	if (0 == uSeqCount)
		Quit("Input file '%s' has no sequences", g_pstrInFileName);
	if (1 == uSeqCount)
		{
		TextFile fileOut(g_pstrOutFileName, true);
		v.ToFile(fileOut);
		return;
		}

	if (uSeqCount > 1)
		MHackStart(v);

// First iteration
	Tree GuideTree;
	if (0 != g_pstrUseTreeFileName)
		{
	// Discourage users...
		if (!g_bUseTreeNoWarn)
			fprintf(stderr, "%s", g_strUseTreeWarning);

	// Read tree from file
		TextFile TreeFile(g_pstrUseTreeFileName);
		GuideTree.FromFile(TreeFile);

	// Make sure tree is rooted
		if (!GuideTree.IsRooted())
			Quit("User tree must be rooted");

		if (GuideTree.GetLeafCount() != uSeqCount)
			Quit("User tree does not match input sequences");

		const unsigned uNodeCount = GuideTree.GetNodeCount();
		for (unsigned uNodeIndex = 0; uNodeIndex < uNodeCount; ++uNodeIndex)
			{
			if (!GuideTree.IsLeaf(uNodeIndex))
				continue;
			const char *LeafName = GuideTree.GetLeafName(uNodeIndex);
			unsigned uSeqIndex;
			bool SeqFound = v.FindName(LeafName, &uSeqIndex);
			if (!SeqFound)
				Quit("Label %s in tree does not match sequences", LeafName);
			unsigned uId = v.GetSeqIdFromName(LeafName);
			GuideTree.SetLeafId(uNodeIndex, uId);
			}
		}
	else
		TreeFromSeqVect(v, GuideTree, g_Cluster1, g_Distance1, g_Root1);

	const char *Tree1 = ValueOpt("Tree1");
	if (0 != Tree1)
		{
		TextFile f(Tree1, true);
		GuideTree.ToFile(f);
		if (g_bClusterOnly)
			return;
		}

	SetMuscleTree(GuideTree);
	ValidateMuscleIds(GuideTree);

	MSA msa;
	ProgressiveAlign(v, GuideTree, msa);
	SetCurrentAlignment(msa);

	if (0 != g_pstrComputeWeightsFileName)
		{
		extern void OutWeights(const char *FileName, const MSA &msa);
		SetMSAWeightsMuscle(msa);
		OutWeights(g_pstrComputeWeightsFileName, msa);
		return;
		}

	ValidateMuscleIds(msa);

	if (1 == g_uMaxIters || 2 == uSeqCount)
		{
		MuscleOutput(msa);
		return;
		}

// Second iteration
	if (0 == g_pstrUseTreeFileName)
		{
		g_bDiags = g_bDiags2;
		SetIter(2);

		RefineTree(msa, GuideTree);

		const char *Tree2 = ValueOpt("Tree2");
		if (0 != Tree2)
			{
			TextFile f(Tree2, true);
			GuideTree.ToFile(f);
			}
		}
////////////////////////////////////////////////////////////////////////////////
//  Experimental third iteration
////////////////////////////////////////////////////////////////////////////////
	{
	//Tree Tree3;
	//TreeFromMSA(msa, Tree3, CLUSTER_UPGMB, DISTANCE_ScoreDist, ROOT_FromClustering);
	//ProgressiveAlign(v, Tree3, msa);
	//MuscleOutput(msa);
	//return;
	}
////////////////////////////////////////////////////////////////////////////////
	SetSeqWeightMethod(g_SeqWeight2);
	SetMuscleTree(GuideTree);

	if (g_bAnchors)
		RefineVert(msa, GuideTree, g_uMaxIters - 2);
	else
		RefineHoriz(msa, GuideTree, g_uMaxIters - 2, false, false);

	ValidateMuscleIds(msa);
	ValidateMuscleIds(GuideTree);

	MuscleOutput(msa);
	Log("Finished %s\n", GetTimeAsStr());
	}
