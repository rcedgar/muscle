#include "muscle.h"
#include "msa.h"
#include "profile.h"
#include "objscore.h"

#define TRACE			0
#define TRACE_SEQPAIR	0
#define TEST_SPFAST		0

extern SCOREMATRIX VTML_LA;
extern SCOREMATRIX PAM200;
extern SCOREMATRIX PAM200NoCenter;
extern SCOREMATRIX VTML_SP;
extern SCOREMATRIX VTML_SPNoCenter;
extern SCOREMATRIX NUC_SP;

SCORE g_SPScoreLetters;
SCORE g_SPScoreGaps;

static SCORE TermGapScore(bool Gap)
	{
	switch (g_TermGaps)
		{
	case TERMGAPS_Full:
		return 0;

	case TERMGAPS_Half:
		if (Gap)
			return g_scoreGapOpen/2;
		return 0;

	case TERMGAPS_Ext:
		if (Gap)
			return g_scoreGapExtend;
		return 0;
		}
	Quit("TermGapScore?!");
	return 0;
	}

SCORE ScoreSeqPairLetters(const MSA &msa1, unsigned uSeqIndex1,
  const MSA &msa2, unsigned uSeqIndex2)
	{
	const unsigned uColCount = msa1.GetColCount();
	const unsigned uColCount2 = msa2.GetColCount();
	if (uColCount != uColCount2)
		Quit("ScoreSeqPairLetters, different lengths");

#if	TRACE_SEQPAIR
	{
	Log("\n");
	Log("ScoreSeqPairLetters\n");
	MSA msaTmp;
	msaTmp.SetSize(2, uColCount);
	msaTmp.CopySeq(0, msa1, uSeqIndex1);
	msaTmp.CopySeq(1, msa2, uSeqIndex2);
	msaTmp.LogMe();
	}
#endif

	SCORE scoreLetters = 0;
	SCORE scoreGaps = 0;
	bool bGapping1 = false;
	bool bGapping2 = false;

	unsigned uColStart = 0;
	bool bLeftTermGap = false;
	for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
		{
		bool bGap1 = msa1.IsGap(uSeqIndex1, uColIndex);
		bool bGap2 = msa2.IsGap(uSeqIndex2, uColIndex);
		if (!bGap1 || !bGap2)
			{
			if (bGap1 || bGap2)
				bLeftTermGap = true;
			uColStart = uColIndex;
			break;
			}
		}

	unsigned uColEnd = uColCount - 1;
	bool bRightTermGap = false;
	for (int iColIndex = (int) uColCount - 1; iColIndex >= 0; --iColIndex)
		{
		bool bGap1 = msa1.IsGap(uSeqIndex1, iColIndex);
		bool bGap2 = msa2.IsGap(uSeqIndex2, iColIndex);
		if (!bGap1 || !bGap2)
			{
			if (bGap1 || bGap2)
				bRightTermGap = true;
			uColEnd = (unsigned) iColIndex;
			break;
			}
		}

#if	TRACE_SEQPAIR
	Log("LeftTermGap=%d RightTermGap=%d\n", bLeftTermGap, bRightTermGap);
#endif

	for (unsigned uColIndex = uColStart; uColIndex <= uColEnd; ++uColIndex)
		{
		unsigned uLetter1 = msa1.GetLetterEx(uSeqIndex1, uColIndex);
		if (uLetter1 >= g_AlphaSize)
			continue;
		unsigned uLetter2 = msa2.GetLetterEx(uSeqIndex2, uColIndex);
		if (uLetter2 >= g_AlphaSize)
			continue;

		SCORE scoreMatch = (*g_ptrScoreMatrix)[uLetter1][uLetter2];
		scoreLetters += scoreMatch;
		}
	return scoreLetters;
	}

SCORE ScoreSeqPairGaps(const MSA &msa1, unsigned uSeqIndex1,
  const MSA &msa2, unsigned uSeqIndex2)
	{
	const unsigned uColCount = msa1.GetColCount();
	const unsigned uColCount2 = msa2.GetColCount();
	if (uColCount != uColCount2)
		Quit("ScoreSeqPairGaps, different lengths");

#if	TRACE_SEQPAIR
	{
	Log("\n");
	Log("ScoreSeqPairGaps\n");
	MSA msaTmp;
	msaTmp.SetSize(2, uColCount);
	msaTmp.CopySeq(0, msa1, uSeqIndex1);
	msaTmp.CopySeq(1, msa2, uSeqIndex2);
	msaTmp.LogMe();
	}
#endif

	SCORE scoreGaps = 0;
	bool bGapping1 = false;
	bool bGapping2 = false;

	unsigned uColStart = 0;
	bool bLeftTermGap = false;
	for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
		{
		bool bGap1 = msa1.IsGap(uSeqIndex1, uColIndex);
		bool bGap2 = msa2.IsGap(uSeqIndex2, uColIndex);
		if (!bGap1 || !bGap2)
			{
			if (bGap1 || bGap2)
				bLeftTermGap = true;
			uColStart = uColIndex;
			break;
			}
		}

	unsigned uColEnd = uColCount - 1;
	bool bRightTermGap = false;
	for (int iColIndex = (int) uColCount - 1; iColIndex >= 0; --iColIndex)
		{
		bool bGap1 = msa1.IsGap(uSeqIndex1, iColIndex);
		bool bGap2 = msa2.IsGap(uSeqIndex2, iColIndex);
		if (!bGap1 || !bGap2)
			{
			if (bGap1 || bGap2)
				bRightTermGap = true;
			uColEnd = (unsigned) iColIndex;
			break;
			}
		}

#if	TRACE_SEQPAIR
	Log("LeftTermGap=%d RightTermGap=%d\n", bLeftTermGap, bRightTermGap);
#endif

	for (unsigned uColIndex = uColStart; uColIndex <= uColEnd; ++uColIndex)
		{
		bool bGap1 = msa1.IsGap(uSeqIndex1, uColIndex);
		bool bGap2 = msa2.IsGap(uSeqIndex2, uColIndex);

		if (bGap1 && bGap2)
			continue;

		if (bGap1)
			{
			if (!bGapping1)
				{
#if	TRACE_SEQPAIR
				Log("Gap open seq 1 col %d\n", uColIndex);
#endif
				if (uColIndex == uColStart)
					scoreGaps += TermGapScore(true);
				else
					scoreGaps += g_scoreGapOpen;
				bGapping1 = true;
				}
			else
				scoreGaps += g_scoreGapExtend;
			continue;
			}

		else if (bGap2)
			{
			if (!bGapping2)
				{
#if	TRACE_SEQPAIR
				Log("Gap open seq 2 col %d\n", uColIndex);
#endif
				if (uColIndex == uColStart)
					scoreGaps += TermGapScore(true);
				else
					scoreGaps += g_scoreGapOpen;
				bGapping2 = true;
				}
			else
				scoreGaps += g_scoreGapExtend;
			continue;
			}

		bGapping1 = false;
		bGapping2 = false;
		}

	if (bGapping1 || bGapping2)
		{
		scoreGaps -= g_scoreGapOpen;
		scoreGaps += TermGapScore(true);
		}
	return scoreGaps;
	}

// The usual sum-of-pairs objective score: sum the score
// of the alignment of each pair of sequences.
SCORE ObjScoreSP(const MSA &msa, SCORE MatchScore[])
	{
#if	TRACE
	Log("==================ObjScoreSP==============\n");
	Log("msa=\n");
	msa.LogMe();
#endif
	g_SPScoreLetters = 0;
	g_SPScoreGaps = 0;

	if (0 != MatchScore)
		{
		const unsigned uColCount = msa.GetColCount();
		for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
			MatchScore[uColIndex] = 0;
		}

	const unsigned uSeqCount = msa.GetSeqCount();
	SCORE scoreTotal = 0;
	unsigned uPairCount = 0;
#if	TRACE
	Log("Seq1  Seq2     wt1     wt2    Letters         Gaps  Unwt.Score    Wt.Score       Total\n");
	Log("----  ----  ------  ------  ----------  ----------  ----------  ----------  ----------\n");
#endif
	for (unsigned uSeqIndex1 = 0; uSeqIndex1 < uSeqCount; ++uSeqIndex1)
		{
		const WEIGHT w1 = msa.GetSeqWeight(uSeqIndex1);
		for (unsigned uSeqIndex2 = uSeqIndex1 + 1; uSeqIndex2 < uSeqCount; ++uSeqIndex2)
			{
			const WEIGHT w2 = msa.GetSeqWeight(uSeqIndex2);
			const WEIGHT w = w1*w2;

			SCORE scoreLetters = ScoreSeqPairLetters(msa, uSeqIndex1, msa, uSeqIndex2);
			SCORE scoreGaps = ScoreSeqPairGaps(msa, uSeqIndex1, msa, uSeqIndex2);
			SCORE scorePair = scoreLetters + scoreGaps;
			++uPairCount;

			scoreTotal += w*scorePair;

			g_SPScoreLetters += w*scoreLetters;
			g_SPScoreGaps += w*scoreGaps;
#if	TRACE
			Log("%4d  %4d  %6.3f  %6.3f  %10.2f  %10.2f  %10.2f  %10.2f  %10.2f >%s >%s\n",
			  uSeqIndex1,
			  uSeqIndex2,
			  w1,
			  w2,
			  scoreLetters,
			  scoreGaps,
			  scorePair,
			  scorePair*w1*w2,
			  scoreTotal,
			  msa.GetSeqName(uSeqIndex1),
			  msa.GetSeqName(uSeqIndex2));
#endif
			}
		}
#if	TEST_SPFAST
	{
	SCORE f = ObjScoreSPFast(msa);
	Log("Fast  = %.6g\n", f);
	Log("Brute = %.6g\n", scoreTotal);
	if (BTEq(f, scoreTotal))
		Log("Agree\n");
	else
		Log("** DISAGREE **\n");
	}
#endif
//	return scoreTotal / uPairCount;
	return scoreTotal;
	}

// Objective score defined as the dynamic programming score.
// Input is two alignments, which must be of the same length.
// Result is the same profile-profile score that is optimized
// by dynamic programming.
SCORE ObjScoreDP(const MSA &msa1, const MSA &msa2, SCORE MatchScore[])
	{
	const unsigned uColCount = msa1.GetColCount();
	if (msa2.GetColCount() != uColCount)
		Quit("ObjScoreDP, must be same length");

	const unsigned uColCount1 = msa1.GetColCount();
	const unsigned uColCount2 = msa2.GetColCount();

	const ProfPos *PA = ProfileFromMSA(msa1);
	const ProfPos *PB = ProfileFromMSA(msa2);

	return ObjScoreDP_Profs(PA, PB, uColCount1, MatchScore);
	}

SCORE ObjScoreDP_Profs(const ProfPos *PA, const ProfPos *PB, unsigned uColCount,
  SCORE MatchScore[])
	{
//#if	TRACE
//	Log("Profile 1:\n");
//	ListProfile(PA, uColCount, &msa1);
//
//	Log("Profile 2:\n");
//	ListProfile(PB, uColCount, &msa2);
//#endif

	SCORE scoreTotal = 0;

	for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
		{
		const ProfPos &PPA = PA[uColIndex];
		const ProfPos &PPB = PB[uColIndex];

		SCORE scoreGap = 0;
		SCORE scoreMatch = 0;
	// If gapped column...
		if (PPA.m_bAllGaps && PPB.m_bAllGaps)
			scoreGap = 0;
		else if (PPA.m_bAllGaps)
			{
			if (uColCount - 1 == uColIndex || !PA[uColIndex+1].m_bAllGaps)
				scoreGap = PPB.m_scoreGapClose;
			if (0 == uColIndex || !PA[uColIndex-1].m_bAllGaps)
				scoreGap += PPB.m_scoreGapOpen;
			//if (0 == scoreGap)
			//	scoreGap = PPB.m_scoreGapExtend;
			}
		else if (PPB.m_bAllGaps)
			{
			if (uColCount - 1 == uColIndex || !PB[uColIndex+1].m_bAllGaps)
				scoreGap = PPA.m_scoreGapClose;
			if (0 == uColIndex || !PB[uColIndex-1].m_bAllGaps)
				scoreGap += PPA.m_scoreGapOpen;
			//if (0 == scoreGap)
			//	scoreGap = PPA.m_scoreGapExtend;
			}
		else
			scoreMatch = ScoreProfPos2(PPA, PPB);

		if (0 != MatchScore)
			MatchScore[uColIndex] = scoreMatch;

		scoreTotal += scoreMatch + scoreGap;

		extern bool g_bTracePPScore;
		extern MSA *g_ptrPPScoreMSA1;
		extern MSA *g_ptrPPScoreMSA2;
		if (g_bTracePPScore)
			{
			const MSA &msa1 = *g_ptrPPScoreMSA1;
			const MSA &msa2 = *g_ptrPPScoreMSA2;
			const unsigned uSeqCount1 = msa1.GetSeqCount();
			const unsigned uSeqCount2 = msa2.GetSeqCount();

			for (unsigned n = 0; n < uSeqCount1; ++n)
				Log("%c", msa1.GetChar(n, uColIndex));
			Log("  ");
			for (unsigned n = 0; n < uSeqCount2; ++n)
				Log("%c", msa2.GetChar(n, uColIndex));
			Log("  %10.3f", scoreMatch);
			if (scoreGap != 0)
				Log("  %10.3f", scoreGap);
			Log("\n");
			}
		}

	delete[] PA;
	delete[] PB;

	return scoreTotal;
	}

// Objective score defined as the sum of profile-sequence
// scores for each sequence in the alignment. The profile
// is computed from the entire alignment, so this includes
// the score of each sequence against itself. This is to
// avoid recomputing the profile each time, so we reduce
// complexity but introduce a questionable approximation.
// The goal is to see if we can exploit the apparent
// improvement in performance of log-expectation score
// over the usual sum-of-pairs by optimizing this
// objective score in the iterative refinement stage.
SCORE ObjScorePS(const MSA &msa, SCORE MatchScore[])
	{
	if (g_PPScore != PPSCORE_LE)
		Quit("FastScoreMSA_LASimple: LA");

	const unsigned uSeqCount = msa.GetSeqCount();
	const unsigned uColCount = msa.GetColCount();

	const ProfPos *Prof = ProfileFromMSA(msa);

	if (0 != MatchScore)
		for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
			MatchScore[uColIndex] = 0;

	SCORE scoreTotal = 0;
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		const WEIGHT weightSeq = msa.GetSeqWeight(uSeqIndex);
		SCORE scoreSeq = 0;
		for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
			{
			const ProfPos &PP = Prof[uColIndex];
			if (msa.IsGap(uSeqIndex, uColIndex))
				{
				bool bOpen = (0 == uColIndex ||
				  !msa.IsGap(uSeqIndex, uColIndex - 1));
				bool bClose = (uColCount - 1 == uColIndex ||
				  !msa.IsGap(uSeqIndex, uColIndex + 1));

				if (bOpen)
					scoreSeq += PP.m_scoreGapOpen;
				if (bClose)
					scoreSeq += PP.m_scoreGapClose;
				//if (!bOpen && !bClose)
				//	scoreSeq += PP.m_scoreGapExtend;
				}
			else if (msa.IsWildcard(uSeqIndex, uColIndex))
				continue;
			else
				{
				unsigned uLetter = msa.GetLetter(uSeqIndex, uColIndex);
				const SCORE scoreMatch = PP.m_AAScores[uLetter];
				if (0 != MatchScore)
					MatchScore[uColIndex] += weightSeq*scoreMatch;
				scoreSeq += scoreMatch;
				}
			}
		scoreTotal += weightSeq*scoreSeq;
		}

	delete[] Prof;
	return scoreTotal;
	}

// The XP score is the sum of the score of each pair of
// sequences between two profiles which are aligned to
// each other. Notice that for two given profiles aligned
// in different ways, the difference in XP score must be
// the same as the difference in SP score because the
// score of a pair of sequences in one profile doesn't
// depend on the alignment.
SCORE ObjScoreXP(const MSA &msa1, const MSA &msa2)
	{
	const unsigned uColCount1 = msa1.GetColCount();
	const unsigned uColCount2 = msa2.GetColCount();
	if (uColCount1 != uColCount2)
		Quit("ObjScoreXP, alignment lengths differ %u %u", uColCount1, uColCount2);

	const unsigned uSeqCount1 = msa1.GetSeqCount();
	const unsigned uSeqCount2 = msa2.GetSeqCount();

#if	TRACE
	Log("     Score  Weight  Weight       Total\n");
	Log("----------  ------  ------  ----------\n");
#endif

	SCORE scoreTotal = 0;
	unsigned uPairCount = 0;
	for (unsigned uSeqIndex1 = 0; uSeqIndex1 < uSeqCount1; ++uSeqIndex1)
		{
		const WEIGHT w1 = msa1.GetSeqWeight(uSeqIndex1);
		for (unsigned uSeqIndex2 = 0; uSeqIndex2 < uSeqCount2; ++uSeqIndex2)
			{
			const WEIGHT w2 = msa2.GetSeqWeight(uSeqIndex2);
			const WEIGHT w = w1*w2;
			SCORE scoreLetters = ScoreSeqPairLetters(msa1, uSeqIndex1, msa2, uSeqIndex2);
			SCORE scoreGaps = ScoreSeqPairGaps(msa1, uSeqIndex1, msa2, uSeqIndex2);
			SCORE scorePair = scoreLetters + scoreGaps;
			scoreTotal += w1*w2*scorePair;
			++uPairCount;
#if	TRACE
			Log("%10.2f  %6.3f  %6.3f  %10.2f  >%s >%s\n",
			  scorePair,
			  w1,
			  w2,
			  scorePair*w1*w2,
			  msa1.GetSeqName(uSeqIndex1),
			  msa2.GetSeqName(uSeqIndex2));
#endif
			}
		}
	if (0 == uPairCount)
		Quit("0 == uPairCount");

#if	TRACE
	Log("msa1=\n");
	msa1.LogMe();
	Log("msa2=\n");
	msa2.LogMe();
	Log("XP=%g\n", scoreTotal);
#endif
//	return scoreTotal / uPairCount;
	return scoreTotal;
	}
