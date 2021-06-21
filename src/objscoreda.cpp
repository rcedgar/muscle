#include "muscle.h"
#include "msa.h"
#include "profile.h"
#include "objscore.h"

#if	DOUBLE_AFFINE

#define TRACE			0
#define TEST_SPFAST		0

static SCORE GapPenalty(unsigned uLength, bool Term, SCORE g, SCORE e)
	{
	//if (Term)
	//	{
	//	switch (g_TermGap)
	//		{
	//	case TERMGAP_Full:
	//		return g + (uLength - 1)*e;

	//	case TERMGAP_Half:
	//		return g/2 + (uLength - 1)*e;

	//	case TERMGAP_Ext:
	//		return uLength*e;
	//		}
	//	Quit("Bad termgap");
	//	}
	//else
	//	return g + (uLength - 1)*e;
	//return MINUS_INFINITY;
	return g + (uLength - 1)*e;
	}

static SCORE GapPenalty(unsigned uLength, bool Term)
	{
	SCORE s1 = GapPenalty(uLength, Term, g_scoreGapOpen, g_scoreGapExtend);
#if	DOUBLE_AFFINE
	SCORE s2 = GapPenalty(uLength, Term, g_scoreGapOpen2, g_scoreGapExtend2);
	if (s1 > s2)
		return s1;
	return s2;
#else
	return s1;
#endif
	}

static const MSA *g_ptrMSA1;
static const MSA *g_ptrMSA2;
static unsigned g_uSeqIndex1;
static unsigned g_uSeqIndex2;

static void LogGap(unsigned uStart, unsigned uEnd, unsigned uGapLength,
  bool bNTerm, bool bCTerm)
	{
	Log("%16.16s  ", "");
	for (unsigned i = 0; i < uStart; ++i)
		Log(" ");
	unsigned uMyLength = 0;
	for (unsigned i = uStart; i <= uEnd; ++i)
		{
		bool bGap1 = g_ptrMSA1->IsGap(g_uSeqIndex1, i);
		bool bGap2 = g_ptrMSA2->IsGap(g_uSeqIndex2, i);
		if (!bGap1 && !bGap2)
			Quit("Error -- neither gapping");
		if (bGap1 && bGap2)
			Log(".");
		else
			{
			++uMyLength;
			Log("-");
			}
		}
	SCORE s = GapPenalty(uGapLength, bNTerm || bCTerm);
	Log(" L=%d N%d C%d s=%.3g", uGapLength, bNTerm, bCTerm, s);
	Log("\n");
	if (uMyLength != uGapLength)
		Quit("Lengths differ");

	}

static SCORE ScoreSeqPair(const MSA &msa1, unsigned uSeqIndex1,
  const MSA &msa2, unsigned uSeqIndex2, SCORE *ptrLetters, SCORE *ptrGaps)
	{
	g_ptrMSA1 = &msa1;
	g_ptrMSA2 = &msa2;
	g_uSeqIndex1 = uSeqIndex1;
	g_uSeqIndex2 = uSeqIndex2;

	const unsigned uColCount = msa1.GetColCount();
	const unsigned uColCount2 = msa2.GetColCount();
	if (uColCount != uColCount2)
		Quit("ScoreSeqPair, different lengths");

#if	TRACE
	Log("ScoreSeqPair\n");
	Log("%16.16s  ", msa1.GetSeqName(uSeqIndex1));
	for (unsigned i = 0; i < uColCount; ++i)
		Log("%c", msa1.GetChar(uSeqIndex1, i));
	Log("\n");
	Log("%16.16s  ", msa2.GetSeqName(uSeqIndex2));
	for (unsigned i = 0; i < uColCount; ++i)
		Log("%c", msa1.GetChar(uSeqIndex2, i));
	Log("\n");
#endif

	SCORE scoreTotal = 0;

// Substitution scores
	unsigned uFirstLetter1 = uInsane;
	unsigned uFirstLetter2 = uInsane;
	unsigned uLastLetter1 = uInsane;
	unsigned uLastLetter2 = uInsane;
	for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
		{
		bool bGap1 = msa1.IsGap(uSeqIndex1, uColIndex);
		bool bGap2 = msa2.IsGap(uSeqIndex2, uColIndex);
		bool bWildcard1 = msa1.IsWildcard(uSeqIndex1, uColIndex);
		bool bWildcard2 = msa2.IsWildcard(uSeqIndex2, uColIndex);

		if (!bGap1)
			{
			if (uInsane == uFirstLetter1)
				uFirstLetter1 = uColIndex;
			uLastLetter1 = uColIndex;
			}
		if (!bGap2)
			{
			if (uInsane == uFirstLetter2)
				uFirstLetter2 = uColIndex;
			uLastLetter2 = uColIndex;
			}

		if (bGap1 || bGap2 || bWildcard1 || bWildcard2)
			continue;

		unsigned uLetter1 = msa1.GetLetter(uSeqIndex1, uColIndex);
		unsigned uLetter2 = msa2.GetLetter(uSeqIndex2, uColIndex);

		SCORE scoreMatch = (*g_ptrScoreMatrix)[uLetter1][uLetter2];
		scoreTotal += scoreMatch;
#if	TRACE
		Log("%c <-> %c = %7.1f  %10.1f\n",
		  msa1.GetChar(uSeqIndex1, uColIndex),
		  msa2.GetChar(uSeqIndex2, uColIndex),
		  scoreMatch,
		  scoreTotal);
#endif
		}
	
	*ptrLetters = scoreTotal;

// Gap penalties
	unsigned uGapLength = uInsane;
	unsigned uGapStartCol = uInsane;
	bool bGapping1 = false;
	bool bGapping2 = false;

	for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
		{
		bool bGap1 = msa1.IsGap(uSeqIndex1, uColIndex);
		bool bGap2 = msa2.IsGap(uSeqIndex2, uColIndex);

		if (bGap1 && bGap2)
			continue;

		if (bGapping1)
			{
			if (bGap1)
				++uGapLength;
			else
				{
				bGapping1 = false;
				bool bNTerm = (uFirstLetter2 == uGapStartCol);
				bool bCTerm = (uLastLetter2 + 1 == uColIndex);
				SCORE scoreGap = GapPenalty(uGapLength, bNTerm || bCTerm);
				scoreTotal += scoreGap;
#if	TRACE
				LogGap(uGapStartCol, uColIndex - 1, uGapLength, bNTerm, bCTerm);
				Log("GAP         %7.1f  %10.1f\n",
				  scoreGap,
				  scoreTotal);
#endif
				}
			continue;
			}
		else
			{
			if (bGap1)
				{
				uGapStartCol = uColIndex;
				bGapping1 = true;
				uGapLength = 1;
				continue;
				}
			}

		if (bGapping2)
			{
			if (bGap2)
				++uGapLength;
			else
				{
				bGapping2 = false;
				bool bNTerm = (uFirstLetter1 == uGapStartCol);
				bool bCTerm = (uLastLetter1 + 1 == uColIndex);
				SCORE scoreGap = GapPenalty(uGapLength, bNTerm || bCTerm);
				scoreTotal += scoreGap;
#if	TRACE
				LogGap(uGapStartCol, uColIndex - 1, uGapLength, bNTerm, bCTerm);
				Log("GAP         %7.1f  %10.1f\n",
				  scoreGap,
				  scoreTotal);
#endif
				}
			}
		else
			{
			if (bGap2)
				{
				uGapStartCol = uColIndex;
				bGapping2 = true;
				uGapLength = 1;
				}
			}
		}

	if (bGapping1 || bGapping2)
		{
		SCORE scoreGap = GapPenalty(uGapLength, true);
		scoreTotal += scoreGap;
#if	TRACE
		LogGap(uGapStartCol, uColCount - 1, uGapLength, false, true);
		Log("GAP         %7.1f  %10.1f\n",
		  scoreGap,
		  scoreTotal);
#endif
		}
	*ptrGaps = scoreTotal - *ptrLetters;
	return scoreTotal;
	}

// The usual sum-of-pairs objective score: sum the score
// of the alignment of each pair of sequences.
SCORE ObjScoreDA(const MSA &msa, SCORE *ptrLetters, SCORE *ptrGaps)
	{
	const unsigned uSeqCount = msa.GetSeqCount();
	SCORE scoreTotal = 0;
	unsigned uPairCount = 0;
#if	TRACE
	msa.LogMe();
	Log("     Score  Weight  Weight       Total\n");
	Log("----------  ------  ------  ----------\n");
#endif
	SCORE TotalLetters = 0;
	SCORE TotalGaps = 0;
	for (unsigned uSeqIndex1 = 0; uSeqIndex1 < uSeqCount; ++uSeqIndex1)
		{
		const WEIGHT w1 = msa.GetSeqWeight(uSeqIndex1);
		for (unsigned uSeqIndex2 = uSeqIndex1 + 1; uSeqIndex2 < uSeqCount; ++uSeqIndex2)
			{
			const WEIGHT w2 = msa.GetSeqWeight(uSeqIndex2);
			const WEIGHT w = w1*w2;
			SCORE Letters;
			SCORE Gaps;
			SCORE scorePair = ScoreSeqPair(msa, uSeqIndex1, msa, uSeqIndex2,
			  &Letters, &Gaps);
			scoreTotal += w1*w2*scorePair;
			TotalLetters += w1*w2*Letters;
			TotalGaps += w1*w2*Gaps;
			++uPairCount;
#if	TRACE
			Log("%10.2f  %6.3f  %6.3f  %10.2f  %d=%s %d=%s\n",
			  scorePair,
			  w1,
			  w2,
			  scorePair*w1*w2,
			  uSeqIndex1,
			  msa.GetSeqName(uSeqIndex1),
			  uSeqIndex2,
			  msa.GetSeqName(uSeqIndex2));
#endif
			}
		}
	*ptrLetters = TotalLetters;
	*ptrGaps = TotalGaps;
	return scoreTotal;
	}

#endif	// DOUBLE_AFFINE
