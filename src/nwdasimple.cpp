#include "muscle.h"
#include <math.h>
#include "pwpath.h"
#include "profile.h"
#include <stdio.h>

#define	TRACE	0

bool g_bKeepSimpleDP;
SCORE *g_DPM;
SCORE *g_DPD;
SCORE *g_DPE;
SCORE *g_DPI;
SCORE *g_DPJ;
char *g_TBM;
char *g_TBD;
char *g_TBE;
char *g_TBI;
char *g_TBJ;

#if	DOUBLE_AFFINE

static char XlatEdgeType(char c)
	{
	if ('E' == c)
		return 'D';
	if ('J' == c)
		return 'I';
	return c;
	}

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

SCORE NWDASimple(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, PWPath &Path)
	{
	assert(uLengthB > 0 && uLengthA > 0);

	const unsigned uPrefixCountA = uLengthA + 1;
	const unsigned uPrefixCountB = uLengthB + 1;

// Allocate DP matrices
	const size_t LM = uPrefixCountA*uPrefixCountB;
	SCORE *DPL_ = new SCORE[LM];
	SCORE *DPM_ = new SCORE[LM];
	SCORE *DPD_ = new SCORE[LM];
	SCORE *DPE_ = new SCORE[LM];
	SCORE *DPI_ = new SCORE[LM];
	SCORE *DPJ_ = new SCORE[LM];

	char *TBM_ = new char[LM];
	char *TBD_ = new char[LM];
	char *TBE_ = new char[LM];
	char *TBI_ = new char[LM];
	char *TBJ_ = new char[LM];

	memset(TBM_, '?', LM);
	memset(TBD_, '?', LM);
	memset(TBE_, '?', LM);
	memset(TBI_, '?', LM);
	memset(TBJ_, '?', LM);

	DPM(0, 0) = 0;
	DPD(0, 0) = MINUS_INFINITY;
	DPE(0, 0) = MINUS_INFINITY;
	DPI(0, 0) = MINUS_INFINITY;
	DPJ(0, 0) = MINUS_INFINITY;

	DPM(1, 0) = MINUS_INFINITY;
	DPD(1, 0) = PA[0].m_scoreGapOpen;
	DPE(1, 0) = PA[0].m_scoreGapOpen2;
	TBD(1, 0) = 'D';
	TBE(1, 0) = 'E';
	DPI(1, 0) = MINUS_INFINITY;
	DPJ(1, 0) = MINUS_INFINITY;

	DPM(0, 1) = MINUS_INFINITY;
	DPD(0, 1) = MINUS_INFINITY;
	DPE(0, 1) = MINUS_INFINITY;
	DPI(0, 1) = PB[0].m_scoreGapOpen;
	DPJ(0, 1) = PB[0].m_scoreGapOpen2;
	TBI(0, 1) = 'I';
	TBJ(0, 1) = 'J';

// Empty prefix of B is special case
	for (unsigned uPrefixLengthA = 2; uPrefixLengthA < uPrefixCountA; ++uPrefixLengthA)
		{
		DPM(uPrefixLengthA, 0) = MINUS_INFINITY;

		DPD(uPrefixLengthA, 0) = DPD(uPrefixLengthA - 1, 0) + g_scoreGapExtend;
		DPE(uPrefixLengthA, 0) = DPE(uPrefixLengthA - 1, 0) + g_scoreGapExtend2;

		TBD(uPrefixLengthA, 0) = 'D';
		TBE(uPrefixLengthA, 0) = 'E';

		DPI(uPrefixLengthA, 0) = MINUS_INFINITY;
		DPJ(uPrefixLengthA, 0) = MINUS_INFINITY;
		}

// Empty prefix of A is special case
	for (unsigned uPrefixLengthB = 2; uPrefixLengthB < uPrefixCountB; ++uPrefixLengthB)
		{
		DPM(0, uPrefixLengthB) = MINUS_INFINITY;

		DPD(0, uPrefixLengthB) = MINUS_INFINITY;
		DPE(0, uPrefixLengthB) = MINUS_INFINITY;

		DPI(0, uPrefixLengthB) = DPI(0, uPrefixLengthB - 1) + g_scoreGapExtend;
		DPJ(0, uPrefixLengthB) = DPJ(0, uPrefixLengthB - 1) + g_scoreGapExtend2;

		TBI(0, uPrefixLengthB) = 'I';
		TBJ(0, uPrefixLengthB) = 'J';
		}

// Special case to agree with NWFast, no D-I transitions so...
	DPD(uLengthA, 0) = MINUS_INFINITY;
	DPE(uLengthA, 0) = MINUS_INFINITY;
//	DPI(0, uLengthB) = MINUS_INFINITY;
//	DPJ(0, uLengthB) = MINUS_INFINITY;

// ============
// Main DP loop
// ============
	SCORE scoreGapCloseB = MINUS_INFINITY;
	SCORE scoreGapClose2B = MINUS_INFINITY;
	for (unsigned uPrefixLengthB = 1; uPrefixLengthB < uPrefixCountB; ++uPrefixLengthB)
		{
		const ProfPos &PPB = PB[uPrefixLengthB - 1];

		SCORE scoreGapCloseA = MINUS_INFINITY;
		SCORE scoreGapClose2A = MINUS_INFINITY;
		for (unsigned uPrefixLengthA = 1; uPrefixLengthA < uPrefixCountA; ++uPrefixLengthA)
			{
			const ProfPos &PPA = PA[uPrefixLengthA - 1];

			{
		// Match M=LetterA+LetterB
			SCORE scoreLL = ScoreProfPos2(PPA, PPB);
			DPL(uPrefixLengthA, uPrefixLengthB) = scoreLL;

			SCORE scoreMM = DPM(uPrefixLengthA-1, uPrefixLengthB-1);
			SCORE scoreDM = DPD(uPrefixLengthA-1, uPrefixLengthB-1) + scoreGapCloseA;
			SCORE scoreEM = DPE(uPrefixLengthA-1, uPrefixLengthB-1) + scoreGapClose2A;
			SCORE scoreIM = DPI(uPrefixLengthA-1, uPrefixLengthB-1) + scoreGapCloseB;
			SCORE scoreJM = DPJ(uPrefixLengthA-1, uPrefixLengthB-1) + scoreGapClose2B;

			SCORE scoreBest;
			if (scoreMM >= scoreDM && scoreMM >= scoreEM && scoreMM >= scoreIM && scoreMM >= scoreJM)
				{
				scoreBest = scoreMM;
				TBM(uPrefixLengthA, uPrefixLengthB) = 'M';
				}
			else if (scoreDM >= scoreMM && scoreDM >= scoreEM && scoreDM >= scoreIM && scoreDM >= scoreJM)
				{
				scoreBest = scoreDM;
				TBM(uPrefixLengthA, uPrefixLengthB) = 'D';
				}
			else if (scoreEM >= scoreMM && scoreEM >= scoreDM && scoreEM >= scoreIM && scoreEM >= scoreJM)
				{
				scoreBest = scoreEM;
				TBM(uPrefixLengthA, uPrefixLengthB) = 'E';
				}
			else if (scoreIM >= scoreMM && scoreIM >= scoreDM && scoreIM >= scoreEM && scoreIM >= scoreJM)
				{
				scoreBest = scoreIM;
				TBM(uPrefixLengthA, uPrefixLengthB) = 'I';
				}
			else
				{
				assert(scoreJM >= scoreMM && scoreJM >= scoreDM && scoreJM >= scoreEM && scoreJM >= scoreIM);
				scoreBest = scoreJM;
				TBM(uPrefixLengthA, uPrefixLengthB) = 'J';
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

			{
		// Delete E=LetterA+GapB
			SCORE scoreME = DPM(uPrefixLengthA-1, uPrefixLengthB) +
			  PA[uPrefixLengthA-1].m_scoreGapOpen2;
			SCORE scoreEE = DPE(uPrefixLengthA-1, uPrefixLengthB) + g_scoreGapExtend2;

			SCORE scoreBest;
			if (scoreME >= scoreEE)
				{
				scoreBest = scoreME;
				TBE(uPrefixLengthA, uPrefixLengthB) = 'M';
				}
			else
				{
				assert(scoreEE >= scoreME);
				scoreBest = scoreEE;
				TBE(uPrefixLengthA, uPrefixLengthB) = 'E';
				}
			DPE(uPrefixLengthA, uPrefixLengthB) = scoreBest;
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

		// Insert J=GapA+LetterB
			{
			SCORE scoreMJ = DPM(uPrefixLengthA, uPrefixLengthB-1) +
			  PB[uPrefixLengthB - 1].m_scoreGapOpen2;
			SCORE scoreJJ = DPJ(uPrefixLengthA, uPrefixLengthB-1) + g_scoreGapExtend2;

			SCORE scoreBest;
			if (scoreMJ >= scoreJJ)
				{
				scoreBest = scoreMJ;
				TBJ(uPrefixLengthA, uPrefixLengthB) = 'M';
				}
			else 
				{
				assert(scoreJJ > scoreMJ);
				scoreBest = scoreJJ;
				TBJ(uPrefixLengthA, uPrefixLengthB) = 'J';
				}
			DPJ(uPrefixLengthA, uPrefixLengthB) = scoreBest;
			}

			scoreGapCloseA = PPA.m_scoreGapClose;
			scoreGapClose2A = PPA.m_scoreGapClose2;
			}
		scoreGapCloseB = PPB.m_scoreGapClose;
		scoreGapClose2B = PPB.m_scoreGapClose2;
		}

#if TRACE
	Log("\n");
	Log("DA Simple DPL:\n");
	ListDP(DPL_, PA, PB, uPrefixCountA, uPrefixCountB);
	Log("\n");
	Log("DA Simple DPM:\n");
	ListDP(DPM_, PA, PB, uPrefixCountA, uPrefixCountB);
	Log("\n");
	Log("DA Simple DPD:\n");
	ListDP(DPD_, PA, PB, uPrefixCountA, uPrefixCountB);
	Log("\n");
	Log("DA Simple DPE:\n");
	ListDP(DPE_, PA, PB, uPrefixCountA, uPrefixCountB);
	Log("\n");
	Log("DA Simple DPI:\n");
	ListDP(DPI_, PA, PB, uPrefixCountA, uPrefixCountB);
	Log("\n");
	Log("DA Simple DPJ:\n");
	ListDP(DPJ_, PA, PB, uPrefixCountA, uPrefixCountB);
	Log("\n");
	Log("DA Simple TBM:\n");
	ListTB(TBM_, PA, PB, uPrefixCountA, uPrefixCountB);
	Log("\n");
	Log("DA Simple TBD:\n");
	ListTB(TBD_, PA, PB, uPrefixCountA, uPrefixCountB);
	Log("\n");
	Log("DA Simple TBE:\n");
	ListTB(TBE_, PA, PB, uPrefixCountA, uPrefixCountB);
	Log("\n");
	Log("DA Simple TBI:\n");
	ListTB(TBI_, PA, PB, uPrefixCountA, uPrefixCountB);
	Log("\n");
	Log("DA Simple TBJ:\n");
	ListTB(TBJ_, PA, PB, uPrefixCountA, uPrefixCountB);
#endif

// Trace-back
// ==========
	Path.Clear();

// Find last edge
	SCORE M = DPM(uLengthA, uLengthB);
	SCORE D = DPD(uLengthA, uLengthB) + PA[uLengthA-1].m_scoreGapClose;
	SCORE E = DPE(uLengthA, uLengthB) + PA[uLengthA-1].m_scoreGapClose2;
	SCORE I = DPI(uLengthA, uLengthB) + PB[uLengthB-1].m_scoreGapClose;
	SCORE J = DPJ(uLengthA, uLengthB) + PB[uLengthB-1].m_scoreGapClose2;
	char cEdgeType = '?';

	SCORE BestScore = M;
	cEdgeType = 'M';
	if (D > BestScore)
		{
		cEdgeType = 'D';
		BestScore = D;
		}
	if (E > BestScore)
		{
		cEdgeType = 'E';
		BestScore = E;
		}
	if (I > BestScore)
		{
		cEdgeType = 'I';
		BestScore = I;
		}
	if (J > BestScore)
		{
		cEdgeType = 'J';
		BestScore = J;
		}

#if	TRACE
	Log("DA Simple: MAB=%.4g DAB=%.4g EAB=%.4g IAB=%.4g JAB=%.4g best=%c\n",
	  M, D, E, I, J, cEdgeType);
#endif

	unsigned PLA = uLengthA;
	unsigned PLB = uLengthB;
	for (;;)
		{
		PWEdge Edge;
		Edge.cType = XlatEdgeType(cEdgeType);
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

		case 'E':
			assert(PLA > 0);
			cEdgeType = TBE(PLA, PLB);
			--PLA;
			break;

		case 'I':
			assert(PLB > 0);
			cEdgeType = TBI(PLA, PLB);
			--PLB;
			break;
		
		case 'J':
			assert(PLB > 0);
			cEdgeType = TBJ(PLA, PLB);
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
		g_DPE = DPE_;
		g_DPI = DPI_;
		g_DPJ = DPJ_;

		g_TBM = TBM_;
		g_TBD = TBD_;
		g_TBE = TBE_;
		g_TBI = TBI_;
		g_TBJ = TBJ_;
		}
	else
		{
		delete[] DPM_;
		delete[] DPD_;
		delete[] DPE_;
		delete[] DPI_;
		delete[] DPJ_;

		delete[] TBM_;
		delete[] TBD_;
		delete[] TBE_;
		delete[] TBI_;
		delete[] TBJ_;
		}

	return BestScore;
	}

#endif // DOUBLE_AFFINE
