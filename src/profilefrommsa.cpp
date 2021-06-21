#include "muscle.h"
#include "msa.h"
#include "profile.h"

#define TRACE	0

static void LogF(FCOUNT f)
	{
	if (f > -0.00001 && f < 0.00001)
		Log("       ");
	else
		Log("  %5.3f", f);
	}

static const char *LocalScoreToStr(SCORE s)
	{
	static char str[16];
	if (s < -1e10 || s > 1e10)
		return "    *";
	sprintf(str, "%5.1f", s);
	return str;
	}

#if	DOUBLE_AFFINE
void ListProfile(const ProfPos *Prof, unsigned uLength, const MSA *ptrMSA)
	{
	Log("  Pos  Occ     LL     LG     GL     GG     Open  Close  Open2  Clos2\n");
	Log("  ---  ---     --     --     --     --     ----  -----  -----  -----\n");
	for (unsigned n = 0; n < uLength; ++n)
		{
		const ProfPos &PP = Prof[n];
		Log("%5u", n);
		LogF(PP.m_fOcc);
		LogF(PP.m_LL);
		LogF(PP.m_LG);
		LogF(PP.m_GL);
		LogF(PP.m_GG);
		Log("  %s", LocalScoreToStr(-PP.m_scoreGapOpen));
		Log("  %s", LocalScoreToStr(-PP.m_scoreGapClose));
		Log("  %s", LocalScoreToStr(-PP.m_scoreGapOpen2));
		Log("  %s", LocalScoreToStr(-PP.m_scoreGapClose2));
		if (0 != ptrMSA)
			{
			const unsigned uSeqCount = ptrMSA->GetSeqCount();
			Log("  ");
			for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
				Log("%c", ptrMSA->GetChar(uSeqIndex, n));
			}
		Log("\n");
		}

	Log("\n");
	Log("  Pos G");
	for (unsigned n = 0; n < g_AlphaSize; ++n)
		Log("     %c", LetterExToChar(n));
	Log("\n");
	Log("  --- -");
	for (unsigned n = 0; n < g_AlphaSize; ++n)
		Log(" -----");
	Log("\n");

	for (unsigned n = 0; n < uLength; ++n)
		{
		const ProfPos &PP = Prof[n];
		Log("%5u", n);
		if (-1 == PP.m_uResidueGroup)
			Log(" -", PP.m_uResidueGroup);
		else
			Log(" %d", PP.m_uResidueGroup);

		for (unsigned uLetter = 0; uLetter < g_AlphaSize; ++uLetter)
			{
			FCOUNT f = PP.m_fcCounts[uLetter];
			if (f == 0.0)
				Log("      ");
			else
				Log(" %5.3f", f);
			}
		if (0 != ptrMSA)
			{
			const unsigned uSeqCount = ptrMSA->GetSeqCount();
			Log("  ");
			for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
				Log("%c", ptrMSA->GetChar(uSeqIndex, n));
			}
		Log("\n");
		}
	}
#endif	// DOUBLE_AFFINE

#if SINGLE_AFFINE
void ListProfile(const ProfPos *Prof, unsigned uLength, const MSA *ptrMSA)
	{
	Log("  Pos  Occ     LL     LG     GL     GG     Open  Close\n");
	Log("  ---  ---     --     --     --     --     ----  -----\n");
	for (unsigned n = 0; n < uLength; ++n)
		{
		const ProfPos &PP = Prof[n];
		Log("%5u", n);
		LogF(PP.m_fOcc);
		LogF(PP.m_LL);
		LogF(PP.m_LG);
		LogF(PP.m_GL);
		LogF(PP.m_GG);
		Log("  %5.1f", -PP.m_scoreGapOpen);
		Log("  %5.1f", -PP.m_scoreGapClose);
		if (0 != ptrMSA)
			{
			const unsigned uSeqCount = ptrMSA->GetSeqCount();
			Log("  ");
			for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
				Log("%c", ptrMSA->GetChar(uSeqIndex, n));
			}
		Log("\n");
		}

	Log("\n");
	Log("  Pos G");
	for (unsigned n = 0; n < g_AlphaSize; ++n)
		Log("     %c", LetterExToChar(n));
	Log("\n");
	Log("  --- -");
	for (unsigned n = 0; n < g_AlphaSize; ++n)
		Log(" -----");
	Log("\n");

	for (unsigned n = 0; n < uLength; ++n)
		{
		const ProfPos &PP = Prof[n];
		Log("%5u", n);
		if (-1 == PP.m_uResidueGroup)
			Log(" -", PP.m_uResidueGroup);
		else
			Log(" %d", PP.m_uResidueGroup);

		for (unsigned uLetter = 0; uLetter < g_AlphaSize; ++uLetter)
			{
			FCOUNT f = PP.m_fcCounts[uLetter];
			if (f == 0.0)
				Log("      ");
			else
				Log(" %5.3f", f);
			}
		if (0 != ptrMSA)
			{
			const unsigned uSeqCount = ptrMSA->GetSeqCount();
			Log("  ");
			for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
				Log("%c", ptrMSA->GetChar(uSeqIndex, n));
			}
		Log("\n");
		}
	}
#endif

void SortCounts(const FCOUNT fcCounts[], unsigned SortOrder[])
	{
	static unsigned InitialSortOrder[MAX_ALPHA] =
		{
		0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19
		};
	memcpy(SortOrder, InitialSortOrder, g_AlphaSize*sizeof(unsigned));

	bool bAny = true;
	while (bAny)
		{
		bAny = false;
		for (unsigned n = 0; n < g_AlphaSize - 1; ++n)
			{
			unsigned i1 = SortOrder[n];
			unsigned i2 = SortOrder[n+1];
			if (fcCounts[i1] < fcCounts[i2])
				{
				SortOrder[n+1] = i1;
				SortOrder[n] = i2;
				bAny = true;
				}
			}
		}
	}

static unsigned AminoGroupFromFCounts(const FCOUNT fcCounts[])
	{
	bool bAny = false;
	unsigned uConsensusResidueGroup = RESIDUE_GROUP_MULTIPLE;
	for (unsigned uLetter = 0; uLetter < 20; ++uLetter)
		{
		if (0 == fcCounts[uLetter])
			continue;
		const unsigned uResidueGroup = ResidueGroup[uLetter];
		if (bAny)
			{
			if (uResidueGroup != uConsensusResidueGroup)
				return RESIDUE_GROUP_MULTIPLE;
			}
		else
			{
			bAny = true;
			uConsensusResidueGroup = uResidueGroup;
			}
		}
	return uConsensusResidueGroup;
	}

static unsigned NucleoGroupFromFCounts(const FCOUNT fcCounts[])
	{
	bool bAny = false;
	unsigned uConsensusResidueGroup = RESIDUE_GROUP_MULTIPLE;
	for (unsigned uLetter = 0; uLetter < 4; ++uLetter)
		{
		if (0 == fcCounts[uLetter])
			continue;
		const unsigned uResidueGroup = uLetter;
		if (bAny)
			{
			if (uResidueGroup != uConsensusResidueGroup)
				return RESIDUE_GROUP_MULTIPLE;
			}
		else
			{
			bAny = true;
			uConsensusResidueGroup = uResidueGroup;
			}
		}
	return uConsensusResidueGroup;
	}

unsigned ResidueGroupFromFCounts(const FCOUNT fcCounts[])
	{
	switch (g_Alpha)
		{
	case ALPHA_Amino:
		return AminoGroupFromFCounts(fcCounts);

	case ALPHA_DNA:
	case ALPHA_RNA:
		return NucleoGroupFromFCounts(fcCounts);
		}
	Quit("ResidueGroupFromFCounts: bad alpha");
	return 0;
	}

ProfPos *ProfileFromMSA(const MSA &a)
	{
	const unsigned uSeqCount = a.GetSeqCount();
	const unsigned uColCount = a.GetColCount();

// Yuck -- cast away const (inconsistent design here).
	SetMSAWeightsMuscle((MSA &) a);

	ProfPos *Pos = new ProfPos[uColCount];

	unsigned uHydrophobicRunLength = 0;
	for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
		{
		ProfPos &PP = Pos[uColIndex];

		PP.m_bAllGaps = a.IsGapColumn(uColIndex);

		FCOUNT fcGapStart;
		FCOUNT fcGapEnd;
		FCOUNT fcGapExtend;
		FCOUNT fOcc;
		a.GetFractionalWeightedCounts(uColIndex, g_bNormalizeCounts, PP.m_fcCounts,
		  &fcGapStart, &fcGapEnd, &fcGapExtend, &fOcc,
		  &PP.m_LL, &PP.m_LG, &PP.m_GL, &PP.m_GG);
		PP.m_fOcc = fOcc;

		SortCounts(PP.m_fcCounts, PP.m_uSortOrder);

		PP.m_uResidueGroup = ResidueGroupFromFCounts(PP.m_fcCounts);

		for (unsigned i = 0; i < g_AlphaSize; ++i)
			{
			SCORE scoreSum = 0;
			for (unsigned j = 0; j < g_AlphaSize; ++j)
				scoreSum += PP.m_fcCounts[j]*(*g_ptrScoreMatrix)[i][j];
			PP.m_AAScores[i] = scoreSum;
			}

		SCORE sStartOcc = (SCORE) (1.0 - fcGapStart);
		SCORE sEndOcc = (SCORE) (1.0 - fcGapEnd);

		PP.m_fcStartOcc = sStartOcc;
		PP.m_fcEndOcc = sEndOcc;

		PP.m_scoreGapOpen = sStartOcc*g_scoreGapOpen/2;
		PP.m_scoreGapClose = sEndOcc*g_scoreGapOpen/2;
#if	DOUBLE_AFFINE
		PP.m_scoreGapOpen2 = sStartOcc*g_scoreGapOpen2/2;
		PP.m_scoreGapClose2 = sEndOcc*g_scoreGapOpen2/2;
#endif
//		PP.m_scoreGapExtend = (SCORE) ((1.0 - fcGapExtend)*scoreGapExtend);

#if	PAF
		if (ALHPA_Amino == g_Alpha && sStartOcc > 0.5)
			{
			extern SCORE PAFactor(const FCOUNT fcCounts[]);
			SCORE paf = PAFactor(PP.m_fcCounts);
			PP.m_scoreGapOpen *= paf;
			PP.m_scoreGapClose *= paf;
			}
#endif
		}

#if	HYDRO
	if (ALPHA_Amino == g_Alpha)
		Hydro(Pos, uColCount);
#endif

#if	TRACE
	{
	Log("ProfileFromMSA\n");
	ListProfile(Pos, uColCount, &a);
	}
#endif
	return Pos;
	}
