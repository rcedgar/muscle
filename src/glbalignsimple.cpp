#include "muscle.h"
#include <math.h>
#include "pwpath.h"
#include "profile.h"
#include <stdio.h>

#define	TRACE	0

#if	1 // SINGLE_AFFINE

extern bool g_bKeepSimpleDP;
extern SCORE *g_DPM;
extern SCORE *g_DPD;
extern SCORE *g_DPI;
extern char *g_TBM;
extern char *g_TBD;
extern char *g_TBI;

static const char *LocalScoreToStr(SCORE s)
	{
	static char str[16];
	if (s < -100000)
		return "     *";
	sprintf(str, "%6.1f", s);
	return str;
	}

static void ListTB(const char *TBM_, const ProfPos *PA, const ProfPos *PB,
  unsigned uPrefixCountA, unsigned uPrefixCountB)
	{
	Log("        ");
	for (unsigned uPrefixLengthB = 0; uPrefixLengthB < uPrefixCountB; ++uPrefixLengthB)
		{
		char c = ' ';
		if (uPrefixLengthB > 0)
			c = ConsensusChar(PB[uPrefixLengthB - 1]);
		Log(" %4u:%c", uPrefixLengthB, c);
		}
	Log("\n");
	for (unsigned uPrefixLengthA = 0; uPrefixLengthA < uPrefixCountA; ++uPrefixLengthA)
		{
		char c = ' ';
		if (uPrefixLengthA > 0)
			c = ConsensusChar(PA[uPrefixLengthA - 1]);
		Log("%4u:%c  ", uPrefixLengthA, c);
		for (unsigned uPrefixLengthB = 0; uPrefixLengthB < uPrefixCountB; ++uPrefixLengthB)
			Log(" %6c", TBM(uPrefixLengthA, uPrefixLengthB));
		Log("\n");
		}
	}

static void ListDP(const SCORE *DPM_, const ProfPos *PA, const ProfPos *PB,
  unsigned uPrefixCountA, unsigned uPrefixCountB)
	{
	Log("        ");
	for (unsigned uPrefixLengthB = 0; uPrefixLengthB < uPrefixCountB; ++uPrefixLengthB)
		{
		char c = ' ';
		if (uPrefixLengthB > 0)
			c = ConsensusChar(PB[uPrefixLengthB - 1]);
		Log(" %4u:%c", uPrefixLengthB, c);
		}
	Log("\n");
	for (unsigned uPrefixLengthA = 0; uPrefixLengthA < uPrefixCountA; ++uPrefixLengthA)
		{
		char c = ' ';
		if (uPrefixLengthA > 0)
			c = ConsensusChar(PA[uPrefixLengthA - 1]);
		Log("%4u:%c  ", uPrefixLengthA, c);
		for (unsigned uPrefixLengthB = 0; uPrefixLengthB < uPrefixCountB; ++uPrefixLengthB)
			Log(" %s", LocalScoreToStr(DPM(uPrefixLengthA, uPrefixLengthB)));
		Log("\n");
		}
	}

SCORE GlobalAlignSimple(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, PWPath &Path)
	{
	assert(uLengthB > 0 && uLengthA > 0);

	SetTermGaps(PA, uLengthA);
	SetTermGaps(PB, uLengthB);

	const unsigned uPrefixCountA = uLengthA + 1;
	const unsigned uPrefixCountB = uLengthB + 1;

// Allocate DP matrices
	const size_t LM = uPrefixCountA*uPrefixCountB;
	SCORE *DPL_ = new SCORE[LM];
	SCORE *DPM_ = new SCORE[LM];
	SCORE *DPD_ = new SCORE[LM];
	SCORE *DPI_ = new SCORE[LM];

	char *TBM_ = new char[LM];
	char *TBD_ = new char[LM];
	char *TBI_ = new char[LM];

	memset(TBM_, '?', LM);
	memset(TBD_, '?', LM);
	memset(TBI_, '?', LM);

	DPM(0, 0) = 0;
	DPD(0, 0) = MINUS_INFINITY;
	DPI(0, 0) = MINUS_INFINITY;

	DPM(1, 0) = MINUS_INFINITY;
	DPD(1, 0) = PA[0].m_scoreGapOpen;
	TBD(1, 0) = 'D';
	DPI(1, 0) = MINUS_INFINITY;

	DPM(0, 1) = MINUS_INFINITY;
	DPD(0, 1) = MINUS_INFINITY;
	DPI(0, 1) = PB[0].m_scoreGapOpen;
	TBI(0, 1) = 'I';

// Empty prefix of B is special case
	for (unsigned uPrefixLengthA = 2; uPrefixLengthA < uPrefixCountA; ++uPrefixLengthA)
		{
	// M=LetterA+LetterB, impossible with empty prefix
		DPM(uPrefixLengthA, 0) = MINUS_INFINITY;

	// D=LetterA+GapB
		DPD(uPrefixLengthA, 0) = DPD(uPrefixLengthA - 1, 0) + g_scoreGapExtend;
		TBD(uPrefixLengthA, 0) = 'D';

	// I=GapA+LetterB, impossible with empty prefix
		DPI(uPrefixLengthA, 0) = MINUS_INFINITY;
		}

// Empty prefix of A is special case
	for (unsigned uPrefixLengthB = 2; uPrefixLengthB < uPrefixCountB; ++uPrefixLengthB)
		{
	// M=LetterA+LetterB, impossible with empty prefix
		DPM(0, uPrefixLengthB) = MINUS_INFINITY;

	// D=LetterA+GapB, impossible with empty prefix
		DPD(0, uPrefixLengthB) = MINUS_INFINITY;

	// I=GapA+LetterB
		DPI(0, uPrefixLengthB) = DPI(0, uPrefixLengthB - 1) + g_scoreGapExtend;
		TBI(0, uPrefixLengthB) = 'I';
		}

// Special case to agree with NWFast, no D-I transitions so...
	DPD(uLengthA, 0) = MINUS_INFINITY;
//	DPI(0, uLengthB) = MINUS_INFINITY;

// ============
// Main DP loop
// ============
	SCORE scoreGapCloseB = MINUS_INFINITY;
	for (unsigned uPrefixLengthB = 1; uPrefixLengthB < uPrefixCountB; ++uPrefixLengthB)
		{
		const ProfPos &PPB = PB[uPrefixLengthB - 1];

		SCORE scoreGapCloseA = MINUS_INFINITY;
		for (unsigned uPrefixLengthA = 1; uPrefixLengthA < uPrefixCountA; ++uPrefixLengthA)
			{
			const ProfPos &PPA = PA[uPrefixLengthA - 1];

			{
		// Match M=LetterA+LetterB
			SCORE scoreLL = ScoreProfPos2(PPA, PPB);
			DPL(uPrefixLengthA, uPrefixLengthB) = scoreLL;

			SCORE scoreMM = DPM(uPrefixLengthA-1, uPrefixLengthB-1);
			SCORE scoreDM = DPD(uPrefixLengthA-1, uPrefixLengthB-1) + scoreGapCloseA;
			SCORE scoreIM = DPI(uPrefixLengthA-1, uPrefixLengthB-1) + scoreGapCloseB;

			SCORE scoreBest;
			if (scoreMM >= scoreDM && scoreMM >= scoreIM)
				{
				scoreBest = scoreMM;
				TBM(uPrefixLengthA, uPrefixLengthB) = 'M';
				}
			else if (scoreDM >= scoreMM && scoreDM >= scoreIM)
				{
				scoreBest = scoreDM;
				TBM(uPrefixLengthA, uPrefixLengthB) = 'D';
				}
			else 
				{
				assert(scoreIM >= scoreMM && scoreIM >= scoreDM);
				scoreBest = scoreIM;
				TBM(uPrefixLengthA, uPrefixLengthB) = 'I';
				}
			DPM(uPrefixLengthA, uPrefixLengthB) = scoreBest + scoreLL;
			}

			{
		// Delete D=LetterA+GapB
			SCORE scoreMD = DPM(uPrefixLengthA-1, uPrefixLengthB) +
			  PA[uPrefixLengthA-1].m_scoreGapOpen;
			SCORE scoreDD = DPD(uPrefixLengthA-1, uPrefixLengthB) + g_scoreGapExtend;

			SCORE scoreBest;
			if (scoreMD >= scoreDD)
				{
				scoreBest = scoreMD;
				TBD(uPrefixLengthA, uPrefixLengthB) = 'M';
				}
			else
				{
				assert(scoreDD >= scoreMD);
				scoreBest = scoreDD;
				TBD(uPrefixLengthA, uPrefixLengthB) = 'D';
				}
			DPD(uPrefixLengthA, uPrefixLengthB) = scoreBest;
			}

		// Insert I=GapA+LetterB
			{
			SCORE scoreMI = DPM(uPrefixLengthA, uPrefixLengthB-1) +
			  PB[uPrefixLengthB - 1].m_scoreGapOpen;
			SCORE scoreII = DPI(uPrefixLengthA, uPrefixLengthB-1) + g_scoreGapExtend;

			SCORE scoreBest;
			if (scoreMI >= scoreII)
				{
				scoreBest = scoreMI;
				TBI(uPrefixLengthA, uPrefixLengthB) = 'M';
				}
			else 
				{
				assert(scoreII > scoreMI);
				scoreBest = scoreII;
				TBI(uPrefixLengthA, uPrefixLengthB) = 'I';
				}
			DPI(uPrefixLengthA, uPrefixLengthB) = scoreBest;
			}

			scoreGapCloseA = PPA.m_scoreGapClose;
			}
		scoreGapCloseB = PPB.m_scoreGapClose;
		}

#if TRACE
	Log("\n");
	Log("Simple DPL:\n");
	ListDP(DPL_, PA, PB, uPrefixCountA, uPrefixCountB);
	Log("\n");
	Log("Simple DPM:\n");
	ListDP(DPM_, PA, PB, uPrefixCountA, uPrefixCountB);
	Log("\n");
	Log("Simple DPD:\n");
	ListDP(DPD_, PA, PB, uPrefixCountA, uPrefixCountB);
	Log("\n");
	Log("Simple DPI:\n");
	ListDP(DPI_, PA, PB, uPrefixCountA, uPrefixCountB);
	Log("\n");
	Log("Simple TBM:\n");
	ListTB(TBM_, PA, PB, uPrefixCountA, uPrefixCountB);
	Log("\n");
	Log("Simple TBD:\n");
	ListTB(TBD_, PA, PB, uPrefixCountA, uPrefixCountB);
	Log("\n");
	Log("Simple TBI:\n");
	ListTB(TBI_, PA, PB, uPrefixCountA, uPrefixCountB);
#endif

// Trace-back
// ==========
	Path.Clear();

// Find last edge
	SCORE M = DPM(uLengthA, uLengthB);
	SCORE D = DPD(uLengthA, uLengthB) + PA[uLengthA-1].m_scoreGapClose;
	SCORE I = DPI(uLengthA, uLengthB) + PB[uLengthB-1].m_scoreGapClose;
	char cEdgeType = '?';

	SCORE BestScore = MINUS_INFINITY;
	if (M >= D && M >= I)
		{
		cEdgeType = 'M';
		BestScore = M;
		}
	else if (D >= M && D >= I)
		{
		cEdgeType = 'D';
		BestScore = D;
		}
	else 
		{
		assert(I >= M && I >= D);
		cEdgeType = 'I';
		BestScore = I;
		}

#if	TRACE
	Log("Simple: MAB=%.4g DAB=%.4g IAB=%.4g best=%c\n", M, D, I, cEdgeType);
#endif

	unsigned PLA = uLengthA;
	unsigned PLB = uLengthB;
	for (;;)
		{
		PWEdge Edge;
		Edge.cType = cEdgeType;
		Edge.uPrefixLengthA = PLA;
		Edge.uPrefixLengthB = PLB;
#if	TRACE
		Log("Prepend %c%d.%d\n", Edge.cType, PLA, PLB);
#endif
		Path.PrependEdge(Edge);

		switch (cEdgeType)
			{
		case 'M':
			assert(PLA > 0);
			assert(PLB > 0);
			cEdgeType = TBM(PLA, PLB);
			--PLA;
			--PLB;
			break;

		case 'D':
			assert(PLA > 0);
			cEdgeType = TBD(PLA, PLB);
			--PLA;
			break;

		case 'I':
			assert(PLB > 0);
			cEdgeType = TBI(PLA, PLB);
			--PLB;
			break;
		
		default:
			Quit("Invalid edge %c", cEdgeType);
			}
		if (0 == PLA && 0 == PLB)
			break;
		}
	Path.Validate();

//	SCORE Score = TraceBack(PA, uLengthA, PB, uLengthB, DPM_, DPD_, DPI_, Path);

#if	TRACE
	SCORE scorePath = FastScorePath2(PA, uLengthA, PB, uLengthB, Path);
	Path.LogMe();
	Log("Score = %s Path = %s\n", LocalScoreToStr(BestScore), LocalScoreToStr(scorePath));
#endif

	if (g_bKeepSimpleDP)
		{
		g_DPM = DPM_;
		g_DPD = DPD_;
		g_DPI = DPI_;

		g_TBM = TBM_;
		g_TBD = TBD_;
		g_TBI = TBI_;
		}
	else
		{
		delete[] DPM_;
		delete[] DPD_;
		delete[] DPI_;

		delete[] TBM_;
		delete[] TBD_;
		delete[] TBI_;
		}

	return BestScore;
	}

#endif // SINLGLE_AFFINE
