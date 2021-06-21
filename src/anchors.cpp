#include "muscle.h"
#include "msa.h"
#include "objscore.h"

#define	TRACE	0

static void WindowSmooth(const SCORE Score[], unsigned uCount, unsigned uWindowLength,
  SCORE SmoothScore[], double dCeil)
	{
#define	Ceil(x)	((SCORE) ((x) > dCeil ? dCeil : (x)))

	if (1 != uWindowLength%2)
		Quit("WindowSmooth=%u must be odd", uWindowLength);

	if (uCount <= uWindowLength)
		{
		for (unsigned i = 0; i < uCount; ++i)
			SmoothScore[i] = 0;
		return;
		}

	const unsigned w2 = uWindowLength/2;
	for (unsigned i = 0; i < w2; ++i)
		{
		SmoothScore[i] = 0;
		SmoothScore[uCount - i - 1] = 0;
		}

	SCORE scoreWindowTotal = 0;
	for (unsigned i = 0; i < uWindowLength; ++i)
		{
		scoreWindowTotal += Ceil(Score[i]);
		}

	for (unsigned i = w2; ; ++i)
		{
		SmoothScore[i] = scoreWindowTotal/uWindowLength;
		if (i == uCount - w2 - 1)
			break;

		scoreWindowTotal -= Ceil(Score[i - w2]);
		scoreWindowTotal += Ceil(Score[i + w2 + 1]);
		}
#undef Ceil
	}

// Find columns that score above the given threshold.
// A range of scores is defined between the average
// and the maximum. The threshold is a fraction 0.0 .. 1.0
// within that range, where 0.0 is the average score
// and 1.0 is the maximum score.
// "Grade" is by analogy with grading on a curve.
static void FindBestColsGrade(const SCORE Score[], unsigned uCount,
  double dThreshold, unsigned BestCols[], unsigned *ptruBestColCount)
	{
	SCORE scoreTotal = 0;
	for (unsigned uIndex = 0; uIndex < uCount; ++uIndex)
		scoreTotal += Score[uIndex];
	const SCORE scoreAvg = scoreTotal / uCount;

	SCORE scoreMax = MINUS_INFINITY;
	for (unsigned uIndex = 0; uIndex < uCount; ++uIndex)
		if (Score[uIndex] > scoreMax)
			scoreMax = Score[uIndex];

	unsigned uBestColCount = 0;
	for (unsigned uIndex = 0; uIndex < uCount; ++uIndex)
		{
		const SCORE s = Score[uIndex];
		const double dHeight = (s - scoreAvg)/(scoreMax - scoreAvg);
		if (dHeight >= dThreshold)
			{
			BestCols[uBestColCount] = uIndex;
			++uBestColCount;
			}
		}
	*ptruBestColCount = uBestColCount;
	}

// Best col only if all following criteria satisfied:
// (1) Score >= min
// (2) Smoothed score >= min
// (3) No gaps.
static void FindBestColsCombo(const MSA &msa, const SCORE Score[],
  const SCORE SmoothScore[], double dMinScore, double dMinSmoothScore,
  unsigned BestCols[], unsigned *ptruBestColCount)
	{
	const unsigned uColCount = msa.GetColCount();

	unsigned uBestColCount = 0;
	for (unsigned uIndex = 0; uIndex < uColCount; ++uIndex)
		{
		if (Score[uIndex] < dMinScore)
			continue;
		if (SmoothScore[uIndex] < dMinSmoothScore)
			continue;
		if (msa.ColumnHasGap(uIndex))
			continue;
		BestCols[uBestColCount] = uIndex;
		++uBestColCount;
		}
	*ptruBestColCount = uBestColCount;
	}

static void ListBestCols(const MSA &msa, const SCORE Score[], const SCORE SmoothScore[],
  unsigned BestCols[], unsigned uBestColCount)
	{
	const unsigned uColCount = msa.GetColCount();
	const unsigned uSeqCount = msa.GetSeqCount();

	Log("Col  ");
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		Log("%u", uSeqIndex%10);
	Log("  ");

	for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
		{
		Log("%3u  ", uColIndex);
		for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
			Log("%c", msa.GetChar(uSeqIndex, uColIndex));

		Log("  %10.3f", Score[uColIndex]);
		Log("  %10.3f", SmoothScore[uColIndex]);

		for (unsigned i = 0; i < uBestColCount; ++i)
			if (BestCols[i] == uColIndex)
				Log(" <-- Best");
		Log("\n");
		}
	}

// If two best columns are found within a window, choose
// the highest-scoring. If more than two, choose the one
// closest to the center of the window.
static void MergeBestCols(const SCORE Scores[], const unsigned BestCols[],
  unsigned uBestColCount, unsigned uWindowLength, unsigned AnchorCols[],
  unsigned *ptruAnchorColCount)
	{
	unsigned uAnchorColCount = 0;
	for (unsigned n = 0; n < uBestColCount; /* update inside loop */)
		{
		unsigned uBestColIndex = BestCols[n];
		unsigned uCountWithinWindow = 0;
		for (unsigned i = n + 1; i < uBestColCount; ++i)
			{
			unsigned uBestColIndex2 = BestCols[i];
			if (uBestColIndex2 - uBestColIndex >= uWindowLength)
				break;
			++uCountWithinWindow;
			}
		unsigned uAnchorCol = uBestColIndex;
		if (1 == uCountWithinWindow)
			{
			unsigned uBestColIndex2 = BestCols[n+1];
			if (Scores[uBestColIndex] > Scores[uBestColIndex2])
				uAnchorCol = uBestColIndex;
			else
				uAnchorCol = uBestColIndex2;
			}
		else if (uCountWithinWindow > 1)
			{
			unsigned uWindowCenter = uBestColIndex + uWindowLength/2;
			int iClosestDist = uWindowLength;
			unsigned uClosestCol = uBestColIndex;
			for (unsigned i = n + 1; i < n + uCountWithinWindow; ++i)
				{
				unsigned uColIndex = BestCols[i];
				int iDist = uColIndex - uBestColIndex;
				if (iDist < 0)
					iDist = -iDist;
				if (iDist < iClosestDist)
					{
					uClosestCol = uColIndex;
					iClosestDist = iDist;
					}
				}
			uAnchorCol = uClosestCol;
			}
		AnchorCols[uAnchorColCount] = uAnchorCol;
		++uAnchorColCount;
		n += uCountWithinWindow + 1;
		}
	*ptruAnchorColCount = uAnchorColCount;
	}

void FindAnchorCols(const MSA &msa, unsigned AnchorCols[],
  unsigned *ptruAnchorColCount)
	{
	const unsigned uColCount = msa.GetColCount();
	if (uColCount < 16)
		{
		*ptruAnchorColCount = 0;
		return;
		}

	SCORE *MatchScore = new SCORE[uColCount];
	SCORE *SmoothScore = new SCORE[uColCount];
	unsigned *BestCols = new unsigned[uColCount];

	GetLetterScores(msa, MatchScore);
	WindowSmooth(MatchScore, uColCount, g_uSmoothWindowLength, SmoothScore,
	  g_dSmoothScoreCeil);

	unsigned uBestColCount;
	FindBestColsCombo(msa, MatchScore, SmoothScore, g_dMinBestColScore, g_dMinSmoothScore,
	  BestCols, &uBestColCount);

#if	TRACE
	ListBestCols(msa, MatchScore, SmoothScore, BestCols, uBestColCount);
#endif

	MergeBestCols(MatchScore, BestCols, uBestColCount, g_uAnchorSpacing, AnchorCols,
	  ptruAnchorColCount);

	delete[] MatchScore;
	delete[] SmoothScore;
	delete[] BestCols;
	}
