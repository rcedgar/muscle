#include "muscle.h"
#include "textfile.h"
#include "msa.h"
#include "tree.h"
#include "profile.h"
#include "objscore.h"

bool g_bTracePPScore = false;
MSA *g_ptrPPScoreMSA1 = 0;
MSA *g_ptrPPScoreMSA2 = 0;

static ProfPos *ProfileFromMSALocal(MSA &msa, Tree &tree)
	{
	const unsigned uSeqCount = msa.GetSeqCount();
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		msa.SetSeqId(uSeqIndex, uSeqIndex);

	TreeFromMSA(msa, tree, g_Cluster2, g_Distance2, g_Root1);
	SetMuscleTree(tree);
	return ProfileFromMSA(msa);
	}

void PPScore()
	{
	if (0 == g_pstrFileName1 || 0 == g_pstrFileName2)
		Quit("-ppscore needs -in1 and -in2");

	SetSeqWeightMethod(g_SeqWeight1);

	TextFile file1(g_pstrFileName1);
	TextFile file2(g_pstrFileName2);

	MSA msa1;
	MSA msa2;

	msa1.FromFile(file1);
	msa2.FromFile(file2);

	const unsigned uLength1 = msa1.GetColCount();
	const unsigned uLength2 = msa2.GetColCount();

	if (uLength1 != uLength2)
		Quit("Profiles must have the same length");

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

	Tree tree1;
	Tree tree2;
	ProfPos *Prof1 = ProfileFromMSALocal(msa1, tree1);
	ProfPos *Prof2 = ProfileFromMSALocal(msa2, tree2);

	g_bTracePPScore = true;
	g_ptrPPScoreMSA1 = &msa1;
	g_ptrPPScoreMSA2 = &msa2;

	SCORE Score = ObjScoreDP_Profs(Prof1, Prof2, uLength1);

	Log("Score=%.4g\n", Score);
	printf("Score=%.4g\n", Score);
	}
