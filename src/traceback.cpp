#include "muscle.h"
#include "profile.h"
#include "pwpath.h"
#include <math.h>

#define TRACE	0

#define EQ(a, b)	(fabs(a-b) < 0.1)

SCORE TraceBack(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, const SCORE *DPM_, const SCORE *DPD_, const SCORE *DPI_,
  PWPath &Path)
	{
#if	TRACE
	Log("\n");
	Log("TraceBack LengthA=%u LengthB=%u\n", uLengthA, uLengthB);
#endif
	assert(uLengthB > 0 && uLengthA > 0);

	const unsigned uPrefixCountA = uLengthA + 1;
	const unsigned uPrefixCountB = uLengthB + 1;

	Path.Clear();

	unsigned uPrefixLengthA = uLengthA;
	unsigned uPrefixLengthB = uLengthB;

	const SCORE scoreM = DPM(uPrefixLengthA, uPrefixLengthB);
	SCORE scoreD = DPD(uPrefixLengthA, uPrefixLengthB);
	SCORE scoreI = DPI(uPrefixLengthA, uPrefixLengthB);

	const ProfPos &LastPPA = PA[uLengthA - 1];
	const ProfPos &LastPPB = PB[uLengthB - 1];

	scoreD += LastPPA.m_scoreGapClose;
	scoreI += LastPPB.m_scoreGapClose;

	char cEdgeType = cInsane;
	SCORE scoreMax;
	if (scoreM >= scoreD && scoreM >= scoreI)
		{
		scoreMax = scoreM;
		cEdgeType = 'M';
		}
	else if (scoreD >= scoreM && scoreD >= scoreI)
		{
		scoreMax = scoreD;
		cEdgeType = 'D';
		}
	else
		{
		assert(scoreI >= scoreM && scoreI >= scoreD);
		scoreMax = scoreI;
		cEdgeType = 'I';
		}

	for (;;)
		{
		if ('S' == cEdgeType)
			break;

		PWEdge Edge;
		Edge.cType = cEdgeType;
		Edge.uPrefixLengthA = uPrefixLengthA;
		Edge.uPrefixLengthB = uPrefixLengthB;
		Path.PrependEdge(Edge);

		char cPrevEdgeType;
		unsigned uPrevPrefixLengthA = uPrefixLengthA;
		unsigned uPrevPrefixLengthB = uPrefixLengthB;

		switch (cEdgeType)
			{
		case 'M':
			{
			assert(uPrefixLengthA > 0);
			assert(uPrefixLengthB > 0);
			const ProfPos &PPA = PA[uPrefixLengthA - 1];
			const ProfPos &PPB = PB[uPrefixLengthB - 1];

			const SCORE Score = DPM(uPrefixLengthA, uPrefixLengthB);
			const SCORE scoreMatch = ScoreProfPos2(PPA, PPB);

			SCORE scoreSM;
			if (1 == uPrefixLengthA && 1 == uPrefixLengthB)
				scoreSM = scoreMatch;
			else
				scoreSM = MINUS_INFINITY;

			SCORE scoreMM = MINUS_INFINITY;
			SCORE scoreDM = MINUS_INFINITY;
			SCORE scoreIM = MINUS_INFINITY;
			if (uPrefixLengthA > 1 && uPrefixLengthB > 1)
				scoreMM = DPM(uPrefixLengthA-1, uPrefixLengthB-1) + scoreMatch;
			if (uPrefixLengthA > 1)
				{
				SCORE scoreTransDM = PA[uPrefixLengthA-2].m_scoreGapClose;
				scoreDM = DPD(uPrefixLengthA-1, uPrefixLengthB-1) + scoreTransDM + scoreMatch;
				}
			if (uPrefixLengthB > 1)
				{
				SCORE scoreTransIM = PB[uPrefixLengthB-2].m_scoreGapClose;
				scoreIM = DPI(uPrefixLengthA-1, uPrefixLengthB-1) + scoreTransIM + scoreMatch;
				}

			if (EQ(scoreMM, Score))
				cPrevEdgeType = 'M';
			else if (EQ(scoreDM, Score))
				cPrevEdgeType = 'D';
			else if (EQ(scoreIM, Score))
				cPrevEdgeType = 'I';
			else if (EQ(scoreSM, Score))
				cPrevEdgeType = 'S';
			else
				Quit("TraceBack: failed to match M score=%g M=%g D=%g I=%g S=%g",
				  Score, scoreMM, scoreDM, scoreIM, scoreSM);

			--uPrevPrefixLengthA;
			--uPrevPrefixLengthB;
			break;
			}

		case 'D':
			{
			assert(uPrefixLengthA > 0);
			const SCORE Score = DPD(uPrefixLengthA, uPrefixLengthB);

			SCORE scoreMD = MINUS_INFINITY;
			SCORE scoreDD = MINUS_INFINITY;
			SCORE scoreSD = MINUS_INFINITY;
			if (uPrefixLengthB == 0)
				{
				if (uPrefixLengthA == 1)
					scoreSD = PA[0].m_scoreGapOpen;
				else
					scoreSD = DPD(uPrefixLengthA - 1, 0);
				}
			if (uPrefixLengthA > 1)
				{
				const ProfPos &PPA = PA[uPrefixLengthA - 1];
				SCORE scoreTransMD = PPA.m_scoreGapOpen;
				scoreMD = DPM(uPrefixLengthA-1, uPrefixLengthB) + scoreTransMD;
				scoreDD = DPD(uPrefixLengthA-1, uPrefixLengthB);
				}

			if (EQ(Score, scoreMD))
				cPrevEdgeType = 'M';
			else if (EQ(Score, scoreDD))
				cPrevEdgeType = 'D';
			else if (EQ(Score, scoreSD))
				cPrevEdgeType = 'S';
			else
				Quit("TraceBack: failed to match D");

			--uPrevPrefixLengthA;
			break;
			}

		case 'I':
			{
			assert(uPrefixLengthB > 0);
			const SCORE Score = DPI(uPrefixLengthA, uPrefixLengthB);

			SCORE scoreMI = MINUS_INFINITY;
			SCORE scoreII = MINUS_INFINITY;
			SCORE scoreSI = MINUS_INFINITY;
			if (uPrefixLengthA == 0)
				{
				if (uPrefixLengthB == 1)
					scoreSI = PB[0].m_scoreGapOpen;
				else
					scoreSI = DPI(0, uPrefixLengthB - 1);
				}
			if (uPrefixLengthB > 1)
				{
				const ProfPos &PPB = PB[uPrefixLengthB - 1];
				SCORE scoreTransMI = PPB.m_scoreGapOpen;
				scoreMI = DPM(uPrefixLengthA, uPrefixLengthB-1) + scoreTransMI;
				scoreII = DPI(uPrefixLengthA, uPrefixLengthB-1);
				}

			if (EQ(Score, scoreMI))
				cPrevEdgeType = 'M';
			else if (EQ(Score, scoreII))
				cPrevEdgeType = 'I';
			else if (EQ(Score, scoreSI))
				cPrevEdgeType = 'S';
			else
				Quit("TraceBack: failed to match I");

			--uPrevPrefixLengthB;
			break;
			}

		default:
			assert(false);
			}
#if	TRACE
		Log("Edge %c%c%u.%u", cPrevEdgeType, cEdgeType, uPrefixLengthA, uPrefixLengthB);
		Log("\n");
#endif
		cEdgeType = cPrevEdgeType;
		uPrefixLengthA = uPrevPrefixLengthA;
		uPrefixLengthB = uPrevPrefixLengthB;
		}

	return scoreMax;
	}
