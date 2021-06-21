#include "muscle.h"
#include "pwpath.h"
#include "profile.h"

#if	DOUBLE_AFFINE

#define	TRACE	0

extern bool g_bKeepSimpleDP;
extern SCORE *g_DPM;
extern SCORE *g_DPD;
extern SCORE *g_DPE;
extern SCORE *g_DPI;
extern SCORE *g_DPJ;
extern char *g_TBM;
extern char *g_TBD;
extern char *g_TBE;
extern char *g_TBI;
extern char *g_TBJ;

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

static void ListDPM(const SCORE *DPM_, const ProfPos *PA, const ProfPos *PB,
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
			{
			SCORE x = (uPrefixLengthA + uPrefixLengthB)*g_scoreGapExtend;
			SCORE s = DPM(uPrefixLengthA, uPrefixLengthB) - x;
			Log(" %s", LocalScoreToStr(s));
			}
		Log("\n");
		}
	}

extern SCORE ScoreProfPos2(const ProfPos &PP, const ProfPos &PPB);

SCORE NWDASimple2(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, PWPath &Path)
	{
	assert(uLengthB > 0 && uLengthA > 0);

	const unsigned uPrefixCountA = uLengthA + 1;
	const unsigned uPrefixCountB = uLengthB + 1;

// Allocate DP matrices
	const size_t LM = uPrefixCountA*uPrefixCountB;
	SCORE *DPM_ = new SCORE[LM];
	SCORE *DPD_ = new SCORE[LM];
	SCORE *DPE_ = new SCORE[LM];
	SCORE *DPI_ = new SCORE[LM];
	SCORE *DPJ_ = new SCORE[LM];
	SCORE *DPL_ = new SCORE[LM];

	char *TBM_ = new char[LM];
	char *TBD_ = new char[LM];
	char *TBE_ = new char[LM];
	char *TBI_ = new char[LM];
	char *TBJ_ = new char[LM];

	memset(DPM_, 0, LM*sizeof(SCORE));
	memset(DPD_, 0, LM*sizeof(SCORE));
	memset(DPE_, 0, LM*sizeof(SCORE));
	memset(DPI_, 0, LM*sizeof(SCORE));
	memset(DPJ_, 0, LM*sizeof(SCORE));

//	memset(DPL_, 0, LM*sizeof(SCORE));

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
	DPI(1, 0) = MINUS_INFINITY;
	DPJ(1, 0) = MINUS_INFINITY;

	DPM(0, 1) = MINUS_INFINITY;
	DPD(0, 1) = MINUS_INFINITY;
	DPE(0, 1) = MINUS_INFINITY;
	DPI(0, 1) = PB[0].m_scoreGapOpen;
	DPJ(0, 1) = PB[0].m_scoreGapOpen2;

// Empty prefix of B is special case
	for (unsigned uPrefixLengthA = 2; uPrefixLengthA < uPrefixCountA; ++uPrefixLengthA)
		{
	// M=LetterA+LetterB, impossible with empty prefix
		DPM(uPrefixLengthA, 0) = MINUS_INFINITY;

	// D=LetterA+GapB
		DPD(uPrefixLengthA, 0) = DPD(uPrefixLengthA - 1, 0) + g_scoreGapExtend;
		TBD(uPrefixLengthA, 0) = 'D';

		DPE(uPrefixLengthA, 0) = DPE(uPrefixLengthA - 1, 0) + g_scoreGapExtend2;
		TBE(uPrefixLengthA, 0) = 'E';

	// I=GapA+LetterB, impossible with empty prefix
		DPI(uPrefixLengthA, 0) = MINUS_INFINITY;
		DPJ(uPrefixLengthA, 0) = MINUS_INFINITY;
		}

// Empty prefix of A is special case
	for (unsigned uPrefixLengthB = 2; uPrefixLengthB < uPrefixCountB; ++uPrefixLengthB)
		{
	// M=LetterA+LetterB, impossible with empty prefix
		DPM(0, uPrefixLengthB) = MINUS_INFINITY;

	// D=LetterA+GapB, impossible with empty prefix
		DPD(0, uPrefixLengthB) = MINUS_INFINITY;
		DPE(0, uPrefixLengthB) = MINUS_INFINITY;

	// I=GapA+LetterB
		DPI(0, uPrefixLengthB) = DPI(0, uPrefixLengthB - 1) + g_scoreGapExtend;
		TBI(0, uPrefixLengthB) = 'I';

		DPJ(0, uPrefixLengthB) = DPJ(0, uPrefixLengthB - 1) + g_scoreGapExtend2;
		TBJ(0, uPrefixLengthB) = 'J';
		}

// ============
// Main DP loop
// ============
	for (unsigned uPrefixLengthB = 1; uPrefixLengthB < uPrefixCountB; ++uPrefixLengthB)
		{
		const ProfPos &PPB = PB[uPrefixLengthB - 1];
		SCORE scoreGapCloseB;
		if (uPrefixLengthB == 1)
			scoreGapCloseB = MINUS_INFINITY;
		else
			scoreGapCloseB = PB[uPrefixLengthB-2].m_scoreGapClose;

		SCORE scoreGapClose2B;
		if (uPrefixLengthB == 1)
			scoreGapClose2B = MINUS_INFINITY;
		else
			scoreGapClose2B = PB[uPrefixLengthB-2].m_scoreGapClose2;

		for (unsigned uPrefixLengthA = 1; uPrefixLengthA < uPrefixCountA; ++uPrefixLengthA)
			{
			const ProfPos &PPA = PA[uPrefixLengthA - 1];

			{
		// Match M=LetterA+LetterB
			SCORE scoreLL = ScoreProfPos2(PPA, PPB);
			DPL(uPrefixLengthA, uPrefixLengthB) = scoreLL;

			SCORE scoreGapCloseA;
			if (uPrefixLengthA == 1)
				scoreGapCloseA = MINUS_INFINITY;
			else
				scoreGapCloseA = PA[uPrefixLengthA-2].m_scoreGapClose;

			SCORE scoreGapClose2A;
			if (uPrefixLengthA == 1)
				scoreGapClose2A = MINUS_INFINITY;
			else
				scoreGapClose2A = PA[uPrefixLengthA-2].m_scoreGapClose2;

			SCORE scoreMM = DPM(uPrefixLengthA-1, uPrefixLengthB-1);
			SCORE scoreDM = DPD(uPrefixLengthA-1, uPrefixLengthB-1) + scoreGapCloseA;
			SCORE scoreEM = DPE(uPrefixLengthA-1, uPrefixLengthB-1) + scoreGapClose2A;
			SCORE scoreIM = DPI(uPrefixLengthA-1, uPrefixLengthB-1) + scoreGapCloseB;
			SCORE scoreJM = DPJ(uPrefixLengthA-1, uPrefixLengthB-1) + scoreGapClose2B;
			SCORE scoreBest;
			if (scoreMM >= scoreDM && scoreMM >= scoreIM && scoreMM >= scoreEM && scoreMM >= scoreJM)
				{
				scoreBest = scoreMM;
				TBM(uPrefixLengthA, uPrefixLengthB) = 'M';
				}
			else if (scoreDM >= scoreMM && scoreDM >= scoreIM && scoreDM >= scoreEM && scoreDM >= scoreJM)
				{
				scoreBest = scoreDM;
				TBM(uPrefixLengthA, uPrefixLengthB) = 'D';
				}
			else if (scoreEM >= scoreMM && scoreEM >= scoreIM && scoreEM >= scoreDM && scoreEM >= scoreJM)
				{
				scoreBest = scoreEM;
				TBM(uPrefixLengthA, uPrefixLengthB) = 'E';
				}
			else if (scoreIM >= scoreMM && scoreIM >= scoreDM && scoreIM >= scoreEM && scoreIM >= scoreJM)
				{
				scoreBest = scoreIM;
				TBM(uPrefixLengthA, uPrefixLengthB) = 'I';
				}
			else if (scoreJM >= scoreMM && scoreJM >= scoreDM && scoreJM >= scoreEM && scoreJM >= scoreIM)
				{
				scoreBest = scoreJM;
				TBM(uPrefixLengthA, uPrefixLengthB) = 'J';
				}
			else
				Quit("Max failed (M)");

			DPM(uPrefixLengthA, uPrefixLengthB) = scoreBest + scoreLL;
			}

			{
		// Delete D=LetterA+GapB
			SCORE scoreMD = DPM(uPrefixLengthA-1, uPrefixLengthB) +
			  PA[uPrefixLengthA-1].m_scoreGapOpen;
			SCORE scoreDD = DPD(uPrefixLengthA-1, uPrefixLengthB) +
			  g_scoreGapExtend;

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
			SCORE scoreEE = DPE(uPrefixLengthA-1, uPrefixLengthB) +
			  g_scoreGapExtend2;

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
			  PB[uPrefixLengthB-1].m_scoreGapOpen;
			SCORE scoreII = DPI(uPrefixLengthA, uPrefixLengthB-1) +
			  g_scoreGapExtend;

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
			  PB[uPrefixLengthB-1].m_scoreGapOpen2;
			SCORE scoreJJ = DPJ(uPrefixLengthA, uPrefixLengthB-1) +
			  g_scoreGapExtend2;

			SCORE scoreBest;
			if (scoreMJ > scoreJJ)
				{
				scoreBest = scoreMJ;
				TBJ(uPrefixLengthA, uPrefixLengthB) = 'M';
				}
			else 
				{
				assert(scoreJJ >= scoreMJ);
				scoreBest = scoreJJ;
				TBJ(uPrefixLengthA, uPrefixLengthB) = 'J';
				}
			DPJ(uPrefixLengthA, uPrefixLengthB) = scoreBest;
			}
			}
		}

// Special case: close gaps at end of alignment
	DPD(uLengthA, uLengthB) += PA[uLengthA-1].m_scoreGapClose;
	DPE(uLengthA, uLengthB) += PA[uLengthA-1].m_scoreGapClose2;

	DPI(uLengthA, uLengthB) += PB[uLengthB-1].m_scoreGapClose;
	DPJ(uLengthA, uLengthB) += PB[uLengthB-1].m_scoreGapClose2;

#if TRACE
	Log("DPL:\n");
	ListDP(DPL_, PA, PB, uPrefixCountA, uPrefixCountB);
	Log("DPM:\n");
	ListDP(DPM_, PA, PB, uPrefixCountA, uPrefixCountB);
	Log("DPD:\n");
	ListDP(DPD_, PA, PB, uPrefixCountA, uPrefixCountB);
	Log("DPE:\n");
	ListDP(DPE_, PA, PB, uPrefixCountA, uPrefixCountB);
	Log("DPI:\n");
	ListDP(DPI_, PA, PB, uPrefixCountA, uPrefixCountB);
	Log("DPJ:\n");
	ListDP(DPJ_, PA, PB, uPrefixCountA, uPrefixCountB);
	Log("TBM:\n");
	ListTB(TBM_, PA, PB, uPrefixCountA, uPrefixCountB);
	Log("TBD:\n");
	ListTB(TBD_, PA, PB, uPrefixCountA, uPrefixCountB);
	Log("TBE:\n");
	ListTB(TBE_, PA, PB, uPrefixCountA, uPrefixCountB);
	Log("TBI:\n");
	ListTB(TBI_, PA, PB, uPrefixCountA, uPrefixCountB);
	Log("TBJ:\n");
	ListTB(TBJ_, PA, PB, uPrefixCountA, uPrefixCountB);
#endif

// ==========
// Trace-back
// ==========

	Path.Clear();

// Find last edge
	char cEdgeType = '?';
	SCORE BestScore = MINUS_INFINITY;
	SCORE M = DPM(uLengthA, uLengthB);
	SCORE D = DPD(uLengthA, uLengthB);
	SCORE E = DPE(uLengthA, uLengthB);
	SCORE I = DPI(uLengthA, uLengthB);
	SCORE J = DPJ(uLengthA, uLengthB);

	if (M >= D && M >= E && M >= I && M >= J)
		{
		cEdgeType = 'M';
		BestScore = M;
		}
	else if (D >= M && D >= E && D >= I && D >= J)
		{
		cEdgeType = 'D';
		BestScore = D;
		}
	else if (E >= M && E >= D && E >= I && E >= J)
		{
		cEdgeType = 'E';
		BestScore = E;
		}
	else if (I >= M && I >= D && I >= E && I >= J)
		{
		cEdgeType = 'I';
		BestScore = I;
		}
	else if (J >= M && J >= D && J >= E && J >= I)
		{
		cEdgeType = 'J';
		BestScore = J;
		}
	else
		Quit("Bad max");

	unsigned PLA = uLengthA;
	unsigned PLB = uLengthB;
	unsigned ECount = 0;
	unsigned JCount = 0;
	for (;;)
		{
#if	TRACE
		Log("TraceBack: %c%u.%u\n", cEdgeType, PLA, PLB);
#endif
		PWEdge Edge;
		Edge.cType = XlatEdgeType(cEdgeType);
		Edge.uPrefixLengthA = PLA;
		Edge.uPrefixLengthB = PLB;
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
			++ECount;
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
			++JCount;
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
	//if (ECount > 0 || JCount > 0)
	//	fprintf(stderr, "E=%d J=%d\n", ECount, JCount);
	Path.Validate();
	if (Path.GetMatchCount() + Path.GetDeleteCount() != uLengthA)
		Quit("Path count A");
	if (Path.GetMatchCount() + Path.GetInsertCount() != uLengthB)
		Quit("Path count B");

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

#if	TRACE
	Log("BestScore=%.6g\n", BestScore);
#endif
	return BestScore;
	}

#endif	// DOUBLE_AFFINE
