#include "muscle.h"
#include "profile.h"
#include "pwpath.h"

SCORE FastScorePath2(const ProfPos *PA, unsigned uLengthA,
  const ProfPos *PB, unsigned uLengthB, const PWPath &Path)
	{
	const unsigned uEdgeCount = Path.GetEdgeCount();
	Log("Edge  SS     PLA   PLB   Match     Gap    Total\n");
	Log("----  --     ---   ---   -----     ---    -----\n");
	char cType = 'S';
	SCORE scoreTotal = 0;
	for (unsigned uEdgeIndex = 0; uEdgeIndex < uEdgeCount; ++uEdgeIndex)
		{
		const PWEdge &Edge = Path.GetEdge(uEdgeIndex);
		const char cPrevType = cType;
		cType = Edge.cType;
		const unsigned uPrefixLengthA = Edge.uPrefixLengthA;
		const unsigned uPrefixLengthB = Edge.uPrefixLengthB;
		bool bGap = false;
		bool bMatch = false;
		SCORE scoreGap = 0;
		SCORE scoreMatch = 0;

		switch (cType)
			{
		case 'M':
			{
			if (0 == uPrefixLengthA || 0 == uPrefixLengthB)
				Quit("FastScorePath2, M zero length");

			const ProfPos &PPA = PA[uPrefixLengthA - 1];
			const ProfPos &PPB = PB[uPrefixLengthB - 1];

			bMatch = true;
			scoreMatch = ScoreProfPos2(PPA, PPB);

			if ('D' == cPrevType)
				{
				bGap = true;
				assert(uPrefixLengthA > 1);
				scoreGap = PA[uPrefixLengthA-2].m_scoreGapClose;
				}
			else if ('I' == cPrevType)
				{
				bGap = true;
				assert(uPrefixLengthB > 1);
				scoreGap = PB[uPrefixLengthB-2].m_scoreGapClose;
				}
			break;
			}

		case 'D':
			{
			if (0 == uPrefixLengthA)
				Quit("FastScorePath2, D zero length");

			const ProfPos &PPA = PA[uPrefixLengthA - 1];
			bGap = true;
			switch (cPrevType)
				{
			case 'S':
				scoreGap = PPA.m_scoreGapOpen;
				break;
			case 'M':
				scoreGap = PPA.m_scoreGapOpen;
				break;
			case 'D':
//				scoreGap = g_scoreGapExtend;
				scoreGap = 0;
				break;
			case 'I':
				Quit("FastScorePath2 DI");
				}
			break;
			}

		case 'I':
			{
			if (0 == uPrefixLengthB)
				Quit("FastScorePath2, I zero length");

			const ProfPos &PPB = PB[uPrefixLengthB - 1];
			bGap = true;
			switch (cPrevType)
				{
			case 'S':
				scoreGap = PPB.m_scoreGapOpen;
				break;
			case 'M':
				scoreGap = PPB.m_scoreGapOpen;
				break;
			case 'I':
				scoreGap = 0;
//				scoreGap = g_scoreGapExtend;
				break;
			case 'D':
				Quit("FastScorePath2 DI");
				}
			break;
			}

		case 'U':
			{
			Quit("FastScorePath2 U");
			}

		default:
			Quit("FastScorePath2: invalid type %c", cType);
			}

		Log("%4u  %c%c  %4u  %4u  ", uEdgeIndex, cPrevType, cType,
		  uPrefixLengthA, uPrefixLengthB);
		if (bMatch)
			Log("%7.1f  ", scoreMatch);
		else
			Log("         ");
		if (bGap)
			Log("%7.1f  ", scoreGap);
		else
			Log("         ");
		SCORE scoreEdge = scoreMatch + scoreGap;
		scoreTotal += scoreEdge;
		Log("%7.1f  %7.1f", scoreEdge, scoreTotal);
		Log("\n");
		}

	SCORE scoreGap = 0;
//	if (!g_bTermGapsHalf)
		switch (cType)
			{
		case 'M':
			scoreGap = 0;
			break;

		case 'D':
			{
			const ProfPos &LastPPA = PA[uLengthA - 1];
			scoreGap = LastPPA.m_scoreGapClose;
			break;
			}

		case 'I':
			{
			const ProfPos &LastPPB = PB[uLengthB - 1];
			scoreGap = LastPPB.m_scoreGapClose;
			break;
			}

		case 'U':
			Quit("Unaligned regions not supported");

		case 'S':
			break;

		default:
			Quit("Invalid type %c", cType);
			}

	Log("      %cE  %4u  %4u           %7.1f\n", cType, uLengthA, uLengthB, scoreGap);
	scoreTotal += scoreGap;

	Log("Total = %g\n", scoreTotal);
	return scoreTotal;
	}
