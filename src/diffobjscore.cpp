#include "muscle.h"
#include "msa.h"
#include "objscore.h"
#include "profile.h"

#define TRACE				0
#define COMPARE_3_52		0
#define BRUTE_LETTERS		0

static SCORE ScoreColLetters(const MSA &msa, unsigned uColIndex)
	{
	SCOREMATRIX &Mx = *g_ptrScoreMatrix;
	const unsigned uSeqCount = msa.GetSeqCount();

#if	BRUTE_LETTERS
	SCORE BruteScore = 0;
	for (unsigned uSeqIndex1 = 0; uSeqIndex1 < uSeqCount; ++uSeqIndex1)
		{
		unsigned uLetter1 = msa.GetLetterEx(uSeqIndex1, uColIndex);
		if (uLetter1 >= g_AlphaSize)
			continue;
		WEIGHT w1 = msa.GetSeqWeight(uSeqIndex1);
		for (unsigned uSeqIndex2 = uSeqIndex1 + 1; uSeqIndex2 < uSeqCount; ++uSeqIndex2)
			{
			unsigned uLetter2 = msa.GetLetterEx(uSeqIndex2, uColIndex);
			if (uLetter2 >= g_AlphaSize)
				continue;
			WEIGHT w2 = msa.GetSeqWeight(uSeqIndex2);
			BruteScore += w1*w2*Mx[uLetter1][uLetter2];
			}
		}
#endif
	
	double N = 0;
	for (unsigned uSeqIndex1 = 0; uSeqIndex1 < uSeqCount; ++uSeqIndex1)
		{
		WEIGHT w = msa.GetSeqWeight(uSeqIndex1);
		N += w;
		}
	if (N <= 0)
		return 0;

	FCOUNT Freqs[20];
	memset(Freqs, 0, sizeof(Freqs));
	SCORE Score = 0;
	for (unsigned uSeqIndex1 = 0; uSeqIndex1 < uSeqCount; ++uSeqIndex1)
		{
		unsigned uLetter = msa.GetLetterEx(uSeqIndex1, uColIndex);
		if (uLetter >= g_AlphaSize)
			continue;
		WEIGHT w = msa.GetSeqWeight(uSeqIndex1);
		Freqs[uLetter] += w;
		Score -= w*w*Mx[uLetter][uLetter];
		}

	for (unsigned uLetter1 = 0; uLetter1 < g_AlphaSize; ++uLetter1)
		{
		const FCOUNT f1 = Freqs[uLetter1];
		Score += f1*f1*Mx[uLetter1][uLetter1];
		for (unsigned uLetter2 = uLetter1 + 1; uLetter2 < g_AlphaSize; ++uLetter2)
			{
			const FCOUNT f2 = Freqs[uLetter2];
			Score += 2*f1*f2*Mx[uLetter1][uLetter2];
			}
		}
	Score /= 2;
#if	BRUTE_LETTERS
	assert(BTEq(BruteScore, Score));
#endif
	return Score;
	}

static SCORE ScoreLetters(const MSA &msa, const unsigned Edges[],
  unsigned uEdgeCount)
	{
	const unsigned uSeqCount = msa.GetSeqCount();
	const unsigned uColCount = msa.GetColCount();

// Letters
	SCORE Score = 0;
	for (unsigned uEdgeIndex = 0; uEdgeIndex < uEdgeCount; ++uEdgeIndex)
		{
		const unsigned uColIndex = Edges[uEdgeIndex];
		assert(uColIndex < uColCount);
		Score += ScoreColLetters(msa, uColIndex);
		}
	return Score;
	}

void GetLetterScores(const MSA &msa, SCORE Scores[])
	{
	const unsigned uColCount = msa.GetColCount();
	for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
		Scores[uColIndex] = ScoreColLetters(msa, uColIndex);
	}

SCORE DiffObjScore(
  const MSA &msa1, const PWPath &Path1, const unsigned Edges1[], unsigned uEdgeCount1, 
  const MSA &msa2, const PWPath &Path2, const unsigned Edges2[], unsigned uEdgeCount2)
	{
#if	TRACE
	{
	Log("============DiffObjScore===========\n");
	Log("msa1:\n");
	msa1.LogMe();
	Log("\n");
	Log("Cols1: ");
	for (unsigned i = 0; i < uEdgeCount1; ++i)
		Log(" %u", Edges1[i]);
	Log("\n\n");
	Log("msa2:\n");
	msa2.LogMe();
	Log("Cols2: ");
	for (unsigned i = 0; i < uEdgeCount2; ++i)
		Log(" %u", Edges2[i]);
	Log("\n\n");
	}
#endif

#if	COMPARE_3_52
	extern SCORE g_SPScoreLetters;
	extern SCORE g_SPScoreGaps;
	SCORE SP1 = ObjScoreSP(msa1);
	SCORE SPLetters1 = g_SPScoreLetters;
	SCORE SPGaps1 = g_SPScoreGaps;

	SCORE SP2 = ObjScoreSP(msa2);
	SCORE SPLetters2 = g_SPScoreLetters;
	SCORE SPGaps2 = g_SPScoreGaps;
	SCORE SPDiffLetters = SPLetters2 - SPLetters1;
	SCORE SPDiffGaps = SPGaps2 - SPGaps1;
	SCORE SPDiff = SPDiffLetters + SPDiffGaps;
#endif

	SCORE Letters1 = ScoreLetters(msa1, Edges1, uEdgeCount1);
	SCORE Letters2 = ScoreLetters(msa2, Edges2, uEdgeCount2);

	SCORE Gaps1 = ScoreGaps(msa1, Edges1, uEdgeCount1);
	SCORE Gaps2 = ScoreGaps(msa2, Edges2, uEdgeCount2);

	SCORE DiffLetters = Letters2 - Letters1;
	SCORE DiffGaps = Gaps2 - Gaps1;
	SCORE Diff = DiffLetters + DiffGaps;

#if	COMPARE_3_52
	Log("ObjScoreSP    Letters1=%.4g  Letters2=%.4g  DiffLetters=%.4g\n",
	  SPLetters1, SPLetters2, SPDiffLetters);

	Log("DiffObjScore  Letters1=%.4g  Letters2=%.4g  DiffLetters=%.4g\n",
	  Letters1, Letters2, DiffLetters);

	Log("ObjScoreSP    Gaps1=%.4g  Gaps2=%.4g  DiffGaps=%.4g\n",
	  SPGaps1, SPGaps2, SPDiffGaps);

	Log("DiffObjScore  Gaps1=%.4g  Gaps2=%.4g  DiffGaps=%.4g\n",
	  Gaps1, Gaps2, DiffGaps);

	Log("SP diff=%.4g DiffObjScore Diff=%.4g\n", SPDiff, Diff);
#endif

	return Diff;
	}
