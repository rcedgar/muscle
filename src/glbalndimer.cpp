#include "muscle.h"
#include <math.h>
#include <stdio.h>	// for sprintf
#include "pwpath.h"
#include "profile.h"
#include "gapscoredimer.h"

#define	TRACE	0

static SCORE TraceBackDimer(  const SCORE *DPM_, const SCORE *DPD_, const SCORE *DPI_,
  const char *TBM_, const char *TBD_, const char *TBI_,
  unsigned uLengthA, unsigned uLengthB, PWPath &Path);

static const char *LocalScoreToStr(SCORE s)
	{
	static char str[16];
	if (MINUS_INFINITY == s)
		return "     *";
	sprintf(str, "%6.3g", s);
	return str;
	}

#if	TRACE
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
		Log("%2d", uPrefixLengthB);
	Log("\n");
	Log("        ");
	for (unsigned uPrefixLengthB = 0; uPrefixLengthB < uPrefixCountB; ++uPrefixLengthB)
		{
		char c = ' ';
		if (uPrefixLengthB > 0)
			c = ConsensusChar(PB[uPrefixLengthB - 1]);
		Log(" %c", c);
		}
	Log("\n");
	for (unsigned uPrefixLengthA = 0; uPrefixLengthA < uPrefixCountA; ++uPrefixLengthA)
		{
		char c = ' ';
		if (uPrefixLengthA > 0)
			c = ConsensusChar(PA[uPrefixLengthA - 1]);
		Log("%4u:%c  ", uPrefixLengthA, c);
		for (unsigned uPrefixLengthB = 0; uPrefixLengthB < uPrefixCountB; ++uPrefixLengthB)
			Log(" %c", TBM(uPrefixLengthA, uPrefixLengthB));
		Log("\n");
		}
	}
#endif // TRACE

static ProfPos PPTerm;
static bool InitializePPTerm()
	{
	PPTerm.m_bAllGaps = false;
	PPTerm.m_LL = 1;
	PPTerm.m_LG = 0;
	PPTerm.m_GL = 0;
	PPTerm.m_GG = 0;
	PPTerm.m_fOcc = 1;
	return true;
	}
static bool PPTermInitialized = InitializePPTerm();

static SCORE ScoreProfPosDimerLE(const ProfPos &PPA, const ProfPos &PPB)
	{
	SCORE Score = 0;
	for (unsigned n = 0; n < 20; ++n)
		{
		const unsigned uLetter = PPA.m_uSortOrder[n];
		const FCOUNT fcLetter = PPA.m_fcCounts[uLetter];
		if (0 == fcLetter)
			break;
		Score += fcLetter*PPB.m_AAScores[uLetter];
		}
	if (0 == Score)
		return -2.5;
	SCORE logScore = logf(Score);
	return (SCORE) (logScore*(PPA.m_fOcc * PPB.m_fOcc));
	}

static SCORE ScoreProfPosDimerPSP(const ProfPos &PPA, const ProfPos &PPB)
	{
	SCORE Score = 0;
	for (unsigned n = 0; n < 20; ++n)
		{
		const unsigned uLetter = PPA.m_uSortOrder[n];
		const FCOUNT fcLetter = PPA.m_fcCounts[uLetter];
		if (0 == fcLetter)
			break;
		Score += fcLetter*PPB.m_AAScores[uLetter];
		}
	return Score;
	}

static SCORE ScoreProfPosDimer(const ProfPos &PPA, const ProfPos &PPB)
	{
	switch (g_PPScore)
		{
	case PPSCORE_LE:
		return ScoreProfPosDimerLE(PPA, PPB);

	case PPSCORE_SP:
	case PPSCORE_SV:
		return ScoreProfPosDimerPSP(PPA, PPB);
		}
	Quit("Invalid g_PPScore");
	return 0;
	}

// Global alignment dynamic programming
// This variant optimizes the profile-profile SP score under the
// dimer approximation.
SCORE GlobalAlignDimer(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, PWPath &Path)
	{
	assert(uLengthB > 0 && uLengthA > 0);

	const unsigned uPrefixCountA = uLengthA + 1;
	const unsigned uPrefixCountB = uLengthB + 1;

// Allocate DP matrices
	const size_t LM = uPrefixCountA*uPrefixCountB;
	SCORE *DPM_ = new SCORE[LM];
	SCORE *DPD_ = new SCORE[LM];
	SCORE *DPI_ = new SCORE[LM];

	char *TBM_ = new char[LM];
	char *TBD_ = new char[LM];
	char *TBI_ = new char[LM];

	DPM(0, 0) = 0;
	DPD(0, 0) = MINUS_INFINITY;
	DPI(0, 0) = MINUS_INFINITY;

	TBM(0, 0) = 'S';
	TBD(0, 0) = '?';
	TBI(0, 0) = '?';

	DPM(1, 0) = MINUS_INFINITY;
	DPD(1, 0) = GapScoreMD(PA[0], PPTerm);
	DPI(1, 0) = MINUS_INFINITY;

	TBM(1, 0) = '?';
	TBD(1, 0) = 'S';
	TBI(1, 0) = '?';

	DPM(0, 1) = MINUS_INFINITY;
	DPD(0, 1) = MINUS_INFINITY;
	DPI(0, 1) = GapScoreMI(PPTerm, PB[0]);

	TBM(0, 1) = '?';
	TBD(0, 1) = '?';
	TBI(0, 1) = 'S';

// Empty prefix of B is special case
	for (unsigned uPrefixLengthA = 2; uPrefixLengthA < uPrefixCountA; ++uPrefixLengthA)
		{
	// M=LetterA+LetterB, impossible with empty prefix
		DPM(uPrefixLengthA, 0) = MINUS_INFINITY;
		TBM(uPrefixLengthA, 0) = '?';

	// D=LetterA+GapB
		DPD(uPrefixLengthA, 0) = DPD(uPrefixLengthA - 1, 0) +
		  GapScoreDD(PA[uPrefixLengthA - 1], PPTerm);
		TBD(uPrefixLengthA, 0) = 'D';

	// I=GapA+LetterB, impossible with empty prefix
		DPI(uPrefixLengthA, 0) = MINUS_INFINITY;
		TBI(uPrefixLengthA, 0) = '?';
		}

// Empty prefix of A is special case
	for (unsigned uPrefixLengthB = 2; uPrefixLengthB < uPrefixCountB; ++uPrefixLengthB)
		{
	// M=LetterA+LetterB, impossible with empty prefix
		DPM(0, uPrefixLengthB) = MINUS_INFINITY;
		TBM(0, uPrefixLengthB) = '?';

	// D=LetterA+GapB, impossible with empty prefix
		DPD(0, uPrefixLengthB) = MINUS_INFINITY;
		TBD(0, uPrefixLengthB) = '?';

	// I=GapA+LetterB
		DPI(0, uPrefixLengthB) = DPI(0, uPrefixLengthB - 1) +
		  GapScoreII(PPTerm, PB[uPrefixLengthB - 1]);
		TBI(0, uPrefixLengthB) = 'I';
		}

// ============
// Main DP loop
// ============
	for (unsigned uPrefixLengthB = 1; uPrefixLengthB < uPrefixCountB; ++uPrefixLengthB)
		{
		const ProfPos &PPB = PB[uPrefixLengthB - 1];
		for (unsigned uPrefixLengthA = 1; uPrefixLengthA < uPrefixCountA; ++uPrefixLengthA)
			{
			const ProfPos &PPA = PA[uPrefixLengthA - 1];
			{
		// Match M=LetterA+LetterB
			SCORE scoreLL = ScoreProfPosDimer(PPA, PPB);

			SCORE scoreMM = DPM(uPrefixLengthA-1, uPrefixLengthB-1) + GapScoreMM(PPA, PPB);
			SCORE scoreDM = DPD(uPrefixLengthA-1, uPrefixLengthB-1) + GapScoreDM(PPA, PPB);
			SCORE scoreIM = DPI(uPrefixLengthA-1, uPrefixLengthB-1) + GapScoreIM(PPA, PPB);

			SCORE scoreBest = scoreMM;
			char c = 'M';
			if (scoreDM > scoreBest)
				{
				scoreBest = scoreDM;
				c = 'D';
				}
			if (scoreIM > scoreBest)
				{
				scoreBest = scoreIM;
				c = 'I';
				}

			DPM(uPrefixLengthA, uPrefixLengthB) = scoreBest + scoreLL;
			TBM(uPrefixLengthA, uPrefixLengthB) = c;
			}
			{
		// Delete D=LetterA+GapB
			SCORE scoreMD = DPM(uPrefixLengthA-1, uPrefixLengthB) + GapScoreMD(PPA, PPB);
			SCORE scoreDD = DPD(uPrefixLengthA-1, uPrefixLengthB) + GapScoreDD(PPA, PPB);
			SCORE scoreID = DPI(uPrefixLengthA-1, uPrefixLengthB) + GapScoreID(PPA, PPB);

			SCORE scoreBest = scoreMD;
			char c = 'M';
			if (scoreDD > scoreBest)
				{
				scoreBest = scoreDD;
				c = 'D';
				}
			if (scoreID > scoreBest)
				{
				scoreBest = scoreID;
				c = 'I';
				}

			DPD(uPrefixLengthA, uPrefixLengthB) = scoreBest;
			TBD(uPrefixLengthA, uPrefixLengthB) = c;
			}
			{
		// Insert I=GapA+LetterB
			SCORE scoreMI = DPM(uPrefixLengthA, uPrefixLengthB-1) + GapScoreMI(PPA, PPB);
			SCORE scoreDI = DPD(uPrefixLengthA, uPrefixLengthB-1) + GapScoreDI(PPA, PPB);
			SCORE scoreII = DPI(uPrefixLengthA, uPrefixLengthB-1) + GapScoreII(PPA, PPB);

			SCORE scoreBest = scoreMI;
			char c = 'M';
			if (scoreDI > scoreBest)
				{
				scoreBest = scoreDI;
				c = 'D';
				}
			if (scoreII > scoreBest)
				{
				scoreBest = scoreII;
				c = 'I';
				}

			DPI(uPrefixLengthA, uPrefixLengthB) = scoreBest;
			TBI(uPrefixLengthA, uPrefixLengthB) = c;
			}
			}
		}

#if TRACE
	Log("DPM:\n");
	ListDP(DPM_, PA, PB, uPrefixCountA, uPrefixCountB);
	Log("DPD:\n");
	ListDP(DPD_, PA, PB, uPrefixCountA, uPrefixCountB);
	Log("DPI:\n");
	ListDP(DPI_, PA, PB, uPrefixCountA, uPrefixCountB);
	Log("TBM:\n");
	ListTB(TBM_, PA, PB, uPrefixCountA, uPrefixCountB);
	Log("TBD:\n");
	ListTB(TBD_, PA, PB, uPrefixCountA, uPrefixCountB);
	Log("TBI:\n");
	ListTB(TBI_, PA, PB, uPrefixCountA, uPrefixCountB);
#endif

	SCORE Score = TraceBackDimer(DPM_, DPD_, DPI_, TBM_, TBD_, TBI_,
	  uLengthA, uLengthB, Path);

#if	TRACE
	Log("GlobalAlignDimer score = %.3g\n", Score);
#endif

	delete[] DPM_;
	delete[] DPD_;
	delete[] DPI_;

	delete[] TBM_;
	delete[] TBD_;
	delete[] TBI_;

	return Score;
	}

static SCORE TraceBackDimer(  const SCORE *DPM_, const SCORE *DPD_, const SCORE *DPI_,
  const char *TBM_, const char *TBD_, const char *TBI_,
  unsigned uLengthA, unsigned uLengthB, PWPath &Path)
	{
	const unsigned uPrefixCountA = uLengthA + 1;

	unsigned uPrefixLengthA = uLengthA;
	unsigned uPrefixLengthB = uLengthB;

	char cEdge = 'M';
	SCORE scoreMax = DPM(uLengthA, uLengthB);
	if (DPD(uLengthA, uLengthB) > scoreMax)
		{
		scoreMax = DPD(uLengthA, uLengthB);
		cEdge = 'D';
		}
	if (DPI(uLengthA, uLengthB) > scoreMax)
		{
		scoreMax = DPI(uLengthA, uLengthB);
		cEdge = 'I';
		}

	for (;;)
		{
		if (0 == uPrefixLengthA && 0 == uPrefixLengthB)
			break;

		PWEdge Edge;
		Edge.cType = cEdge;
		Edge.uPrefixLengthA = uPrefixLengthA;
		Edge.uPrefixLengthB = uPrefixLengthB;
		Path.PrependEdge(Edge);

#if TRACE
		Log("PLA=%u PLB=%u Edge=%c\n", uPrefixLengthA, uPrefixLengthB, cEdge);
#endif
		switch (cEdge)
			{
		case 'M':
			assert(uPrefixLengthA > 0 && uPrefixLengthB > 0);
			cEdge = TBM(uPrefixLengthA, uPrefixLengthB);
			--uPrefixLengthA;
			--uPrefixLengthB;
			break;
		case 'D':
			assert(uPrefixLengthA > 0);
			cEdge = TBD(uPrefixLengthA, uPrefixLengthB);
			--uPrefixLengthA;
			break;
		case 'I':
			assert(uPrefixLengthB > 0);
			cEdge = TBI(uPrefixLengthA, uPrefixLengthB);
			--uPrefixLengthB;
			break;
		default:
			Quit("Invalid edge PLA=%u PLB=%u %c", uPrefixLengthA, uPrefixLengthB, cEdge);
			}
		}
#if	TRACE
	Path.LogMe();
#endif
	return scoreMax;
	}
