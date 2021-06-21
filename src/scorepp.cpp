#include "muscle.h"
#include "profile.h"

char ConsensusChar(const ProfPos &PP)
	{
	unsigned uMostCommonLetter = 0;
	FCOUNT fcMostCommon = PP.m_fcCounts[0];
	bool bMoreThanOneLetter = false;
	bool bAnyLetter = false;
	for (unsigned uLetter = 0; uLetter < g_AlphaSize; ++uLetter)
		{
		const FCOUNT fc = PP.m_fcCounts[uLetter];
		if (fc > 0)
			{
			if (bAnyLetter)
				bMoreThanOneLetter = true;
			bAnyLetter = true;
			}
		if (fc > fcMostCommon)
			{
			uMostCommonLetter = uLetter;
			fcMostCommon = fc;
			}
		}
	if (!bAnyLetter)
		return '-';
	char c = LetterToChar(uMostCommonLetter);
	if (bMoreThanOneLetter)
		return UnalignChar(c);
	return c;
	}

SCORE ScoreProfPos2LA(const ProfPos &PPA, const ProfPos &PPB)
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
	return (SCORE) ((logScore - g_scoreCenter)*(PPA.m_fOcc * PPB.m_fOcc));
	}

SCORE ScoreProfPos2NS(const ProfPos &PPA, const ProfPos &PPB)
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
	return Score - g_scoreCenter;
	}

SCORE ScoreProfPos2SP(const ProfPos &PPA, const ProfPos &PPB)
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
	return Score - g_scoreCenter;
	}

SCORE ScoreProfPos2SPN(const ProfPos &PPA, const ProfPos &PPB)
	{
	SCORE Score = 0;
	for (unsigned n = 0; n < 4; ++n)
		{
		const unsigned uLetter = PPA.m_uSortOrder[n];
		const FCOUNT fcLetter = PPA.m_fcCounts[uLetter];
		if (0 == fcLetter)
			break;
		Score += fcLetter*PPB.m_AAScores[uLetter];
		}
	return Score - g_scoreCenter;
	}

SCORE ScoreProfPos2(const ProfPos &PPA, const ProfPos &PPB)
	{
	if (PPSCORE_SP == g_PPScore)
		return ScoreProfPos2NS(PPA, PPB);
	else if (PPSCORE_LE == g_PPScore)
		return ScoreProfPos2LA(PPA, PPB);
	else if (PPSCORE_SV == g_PPScore)
		return ScoreProfPos2SP(PPA, PPB);
	else if (PPSCORE_SPN == g_PPScore)
		return ScoreProfPos2SPN(PPA, PPB);
	Quit("Invalid g_PPScore");
	return 0;
	}
