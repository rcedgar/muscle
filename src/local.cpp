#include "muscle.h"
#include "textfile.h"
#include "msa.h"
#include "profile.h"
#include "pwpath.h"
#include "tree.h"

#define TRACE	0

static void MSAFromFileName(const char *FileName, MSA &a)
	{
	TextFile File(FileName);
	a.FromFile(File);
	}

static ProfPos *ProfileFromMSALocal(MSA &msa, Tree &tree)
	{
	const unsigned uSeqCount = msa.GetSeqCount();
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		msa.SetSeqId(uSeqIndex, uSeqIndex);

	TreeFromMSA(msa, tree, g_Cluster1, g_Distance1, g_Root1);
	SetMuscleTree(tree);
	return ProfileFromMSA(msa);
	}

void Local()
	{
	if (0 == g_pstrFileName1 || 0 == g_pstrFileName2)
		Quit("Must specify both -in1 and -in2 for -sw");

	SetSeqWeightMethod(g_SeqWeight1);

	MSA msa1;
	MSA msa2;

	MSAFromFileName(g_pstrFileName1, msa1);
	MSAFromFileName(g_pstrFileName2, msa2);

	ALPHA Alpha = ALPHA_Undefined;
	switch (g_SeqType)
		{
	case SEQTYPE_Auto:
		Alpha = msa1.GuessAlpha();
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
		Quit("Invalid SeqType");
		}
	SetAlpha(Alpha);

	msa1.FixAlpha();
	msa2.FixAlpha();

	if (ALPHA_DNA == Alpha || ALPHA_RNA == Alpha)
		SetPPScore(PPSCORE_SPN);

	const unsigned uSeqCount1 = msa1.GetSeqCount();
	const unsigned uSeqCount2 = msa2.GetSeqCount();
	const unsigned uMaxSeqCount = (uSeqCount1 > uSeqCount2 ? uSeqCount1 : uSeqCount2);
	MSA::SetIdCount(uMaxSeqCount);

	unsigned uLength1 = msa1.GetColCount();
	unsigned uLength2 = msa2.GetColCount();

	Tree tree1;
	Tree tree2;

	ProfPos *Prof1 = ProfileFromMSALocal(msa1, tree1);
	ProfPos *Prof2 = ProfileFromMSALocal(msa2, tree2);

	PWPath Path;
	SW(Prof1, uLength1, Prof2, uLength2, Path);

#if	TRACE
	Path.LogMe();
#endif

	MSA msaOut;
	AlignTwoMSAsGivenPathSW(Path, msa1, msa2, msaOut);

#if	TRACE
	msaOut.LogMe();
#endif

	TextFile fileOut(g_pstrOutFileName, true);
	msaOut.ToFile(fileOut);
	}
