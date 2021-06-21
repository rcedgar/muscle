#include "muscle.h"
#include "objscore.h"
#include "msa.h"
#include "textfile.h"
#include "pwpath.h"

const unsigned INDELS = 1;

static void GetPos(const char Str[], unsigned L, int *pi1, int *pi2)
	{
	int i1;
	for (;;)
		{
		i1 = rand()%(L-2) + 1;
		if (Str[i1] == 'M')
			break;
		}
	int i2;
	for (;;)
		{
		i2 = rand()%(L-2) + 1;
		if (i1 != i2 && Str[i2] == 'M')
			break;
		}
	*pi1 = i1;
	*pi2 = i2;
	}

static void MakePath(unsigned uSeqLength, unsigned uIndelCount, char Str[])
	{
	unsigned uPathLength = uSeqLength + uIndelCount;
	for (unsigned i = 0; i < uPathLength; ++i)
		Str[i] = 'M';

	for (unsigned i = 0; i < uIndelCount; ++i)
		{
		int i1, i2;
		GetPos(Str, uPathLength, &i1, &i2);
		Str[i1] = 'D';
		Str[i2] = 'I';
		}

	Str[uPathLength] = 0;
	Log("MakePath=%s\n", Str);
	}

void SPTest()
	{
	SetPPScore(PPSCORE_SV);

	SetListFileName("c:\\tmp\\muscle.log", false);

	TextFile file1("c:\\tmp\\msa1.afa");
	TextFile file2("c:\\tmp\\msa2.afa");

	MSA msa1;
	MSA msa2;

	msa1.FromFile(file1);
	msa2.FromFile(file2);

	Log("msa1=\n");
	msa1.LogMe();
	Log("msa2=\n");
	msa2.LogMe();

	const unsigned uColCount = msa1.GetColCount();
	if (msa2.GetColCount() != uColCount)
		Quit("Different lengths");

	const unsigned uSeqCount1 = msa1.GetSeqCount();
	const unsigned uSeqCount2 = msa2.GetSeqCount();
	const unsigned uSeqCount = uSeqCount1 + uSeqCount2;

	MSA::SetIdCount(uSeqCount);

	for (unsigned uSeqIndex1 = 0; uSeqIndex1 < uSeqCount1; ++uSeqIndex1)
		{
		msa1.SetSeqWeight(uSeqIndex1, 1.0);
		msa1.SetSeqId(uSeqIndex1, uSeqIndex1);
		}

	for (unsigned uSeqIndex2 = 0; uSeqIndex2 < uSeqCount2; ++uSeqIndex2)
		{
		msa2.SetSeqWeight(uSeqIndex2, 1.0);
		msa2.SetSeqId(uSeqIndex2, uSeqCount1 + uSeqIndex2);
		}

	MSA alnA;
	MSA alnB;

	char strPathA[1024];
	char strPathB[1024];
	MakePath(uColCount, INDELS, strPathA);
	MakePath(uColCount, INDELS, strPathB);

	PWPath PathA;
	PWPath PathB;
	PathA.FromStr(strPathA);
	PathB.FromStr(strPathB);

	Log("PathA=\n");
	PathA.LogMe();
	Log("PathB=\n");
	PathB.LogMe();

	AlignTwoMSAsGivenPath(PathA, msa1, msa2, alnA);
	AlignTwoMSAsGivenPath(PathB, msa1, msa2, alnB);

	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		alnA.SetSeqWeight(uSeqIndex, 1.0);
		alnB.SetSeqWeight(uSeqIndex, 1.0);
		}

	unsigned Seqs1[1024];
	unsigned Seqs2[1024];

	for (unsigned uSeqIndex1 = 0; uSeqIndex1 < uSeqCount1; ++uSeqIndex1)
		Seqs1[uSeqIndex1] = uSeqIndex1;

	for (unsigned uSeqIndex2 = 0; uSeqIndex2 < uSeqCount2; ++uSeqIndex2)
		Seqs2[uSeqIndex2] = uSeqCount1 + uSeqIndex2;

	MSA msaA1;
	MSA msaA2;
	MSA msaB1;
	MSA msaB2;
	MSAFromSeqSubset(alnA, Seqs1, uSeqCount1, msaA1);
	MSAFromSeqSubset(alnB, Seqs1, uSeqCount1, msaB1);
	MSAFromSeqSubset(alnA, Seqs2, uSeqCount2, msaA2);
	MSAFromSeqSubset(alnB, Seqs2, uSeqCount2, msaB2);

	for (unsigned uSeqIndex1 = 0; uSeqIndex1 < uSeqCount1; ++uSeqIndex1)
		{
		msaA1.SetSeqWeight(uSeqIndex1, 1.0);
		msaB1.SetSeqWeight(uSeqIndex1, 1.0);
		}

	for (unsigned uSeqIndex2 = 0; uSeqIndex2 < uSeqCount2; ++uSeqIndex2)
		{
		msaA2.SetSeqWeight(uSeqIndex2, 1.0);
		msaB2.SetSeqWeight(uSeqIndex2, 1.0);
		}

	Log("msaA1=\n");
	msaA1.LogMe();

	Log("msaB1=\n");
	msaB1.LogMe();

	Log("msaA2=\n");
	msaA2.LogMe();

	Log("msaB2=\n");
	msaB2.LogMe();

	Log("alnA=\n");
	alnA.LogMe();

	Log("AlnB=\n");
	alnB.LogMe();

	Log("\nSPA\n---\n");
	SCORE SPA = ObjScoreSP(alnA);
	Log("\nSPB\n---\n");
	SCORE SPB = ObjScoreSP(alnB);

	Log("\nXPA\n---\n");
	SCORE XPA = ObjScoreXP(msaA1, msaA2);
	Log("\nXPB\n---\n");
	SCORE XPB = ObjScoreXP(msaB1, msaB2);

	Log("SPA=%.4g SPB=%.4g Diff=%.4g\n", SPA, SPB, SPA - SPB);
	Log("XPA=%.4g XPB=%.4g Diff=%.4g\n", XPA, XPB, XPA - XPB);
	}
