#include "muscle.h"
#include "profile.h"
#include "pwpath.h"

#define	OCC	1

struct DP_MEMORY
	{
	unsigned uLength;
	SCORE *GapOpenA;
	SCORE *GapOpenB;
	SCORE *GapCloseA;
	SCORE *GapCloseB;
	SCORE *MPrev;
	SCORE *MCurr;
	SCORE *MWork;
	SCORE *DPrev;
	SCORE *DCurr;
	SCORE *DWork;
	SCORE **ScoreMxB;
#if	OCC
	FCOUNT *OccA;
	FCOUNT *OccB;
#endif
	unsigned **SortOrderA;
	unsigned *uDeletePos;
	FCOUNT **FreqsA;
	int **TraceBack;
	};

static struct DP_MEMORY DPM;

static void AllocDPMem(unsigned uLengthA, unsigned uLengthB)
	{
// Max prefix length
	unsigned uLength = (uLengthA > uLengthB ? uLengthA : uLengthB) + 1;
	if (uLength < DPM.uLength)
		return;

// Add 256 to allow for future expansion and
// round up to next multiple of 32.
	uLength += 256;
	uLength += 32 - uLength%32;

	const unsigned uOldLength = DPM.uLength;
	if (uOldLength > 0)
		{
		for (unsigned i = 0; i < uOldLength; ++i)
			{
			delete[] DPM.TraceBack[i];
			delete[] DPM.FreqsA[i];
			delete[] DPM.SortOrderA[i];
			}
		for (unsigned n = 0; n < 20; ++n)
			delete[] DPM.ScoreMxB[n];

		delete[] DPM.MPrev;
		delete[] DPM.MCurr;
		delete[] DPM.MWork;
		delete[] DPM.DPrev;
		delete[] DPM.DCurr;
		delete[] DPM.DWork;
		delete[] DPM.uDeletePos;
		delete[] DPM.GapOpenA;
		delete[] DPM.GapOpenB;
		delete[] DPM.GapCloseA;
		delete[] DPM.GapCloseB;
		delete[] DPM.SortOrderA;
		delete[] DPM.FreqsA;
		delete[] DPM.ScoreMxB;
		delete[] DPM.TraceBack;
#if	OCC
		delete[] DPM.OccA;
		delete[] DPM.OccB;
#endif
		}

	DPM.uLength = uLength;

	DPM.GapOpenA = new SCORE[uLength];
	DPM.GapOpenB = new SCORE[uLength];
	DPM.GapCloseA = new SCORE[uLength];
	DPM.GapCloseB = new SCORE[uLength];
#if	OCC
	DPM.OccA = new FCOUNT[uLength];
	DPM.OccB = new FCOUNT[uLength];
#endif

	DPM.SortOrderA = new unsigned*[uLength];
	DPM.FreqsA = new FCOUNT*[uLength];
	DPM.ScoreMxB = new SCORE*[20];
	DPM.MPrev = new SCORE[uLength];
	DPM.MCurr = new SCORE[uLength];
	DPM.MWork = new SCORE[uLength];

	DPM.DPrev = new SCORE[uLength];
	DPM.DCurr = new SCORE[uLength];
	DPM.DWork = new SCORE[uLength];
	DPM.uDeletePos = new unsigned[uLength];

	DPM.TraceBack = new int*[uLength];

	for (unsigned uLetter = 0; uLetter < 20; ++uLetter)
		DPM.ScoreMxB[uLetter] = new SCORE[uLength];

	for (unsigned i = 0; i < uLength; ++i)
		{
		DPM.SortOrderA[i] = new unsigned[20];
		DPM.FreqsA[i] = new FCOUNT[20];
		DPM.TraceBack[i] = new int[uLength];
		}
	}

SCORE GlobalAlignLE(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, PWPath &Path)
	{
	SetTermGaps(PA, uLengthA);
	SetTermGaps(PB, uLengthB);

	const unsigned uPrefixCountA = uLengthA + 1;
	const unsigned uPrefixCountB = uLengthB + 1;

	AllocDPMem(uLengthA, uLengthB);

	SCORE *GapOpenA = DPM.GapOpenA;
	SCORE *GapOpenB = DPM.GapOpenB;
	SCORE *GapCloseA = DPM.GapCloseA;
	SCORE *GapCloseB = DPM.GapCloseB;

	unsigned **SortOrderA = DPM.SortOrderA;
	FCOUNT **FreqsA = DPM.FreqsA;
	SCORE **ScoreMxB = DPM.ScoreMxB;
	SCORE *MPrev = DPM.MPrev;
	SCORE *MCurr = DPM.MCurr;
	SCORE *MWork = DPM.MWork;

	SCORE *DPrev = DPM.DPrev;
	SCORE *DCurr = DPM.DCurr;
	SCORE *DWork = DPM.DWork;

#if	OCC
	FCOUNT *OccA = DPM.OccA;
	FCOUNT *OccB = DPM.OccB;
#endif

	unsigned *uDeletePos = DPM.uDeletePos;

	int **TraceBack = DPM.TraceBack;

	for (unsigned i = 0; i < uLengthA; ++i)
		{
		GapOpenA[i] = PA[i].m_scoreGapOpen;
		GapCloseA[i] = PA[i].m_scoreGapClose;
#if	OCC
		OccA[i] = PA[i].m_fOcc;
#endif

		for (unsigned uLetter = 0; uLetter < 20; ++uLetter)
			{
			SortOrderA[i][uLetter] = PA[i].m_uSortOrder[uLetter];
			FreqsA[i][uLetter] = PA[i].m_fcCounts[uLetter];
			}
		}

	for (unsigned j = 0; j < uLengthB; ++j)
		{
		GapOpenB[j] = PB[j].m_scoreGapOpen;
		GapCloseB[j] = PB[j].m_scoreGapClose;
#if	OCC
		OccB[j] = PB[j].m_fOcc;
#endif
		}

	for (unsigned uLetter = 0; uLetter < 20; ++uLetter)
		{
		for (unsigned j = 0; j < uLengthB; ++j)
			ScoreMxB[uLetter][j] = PB[j].m_AAScores[uLetter];
		}

	for (unsigned i = 0; i < uPrefixCountA; ++i)
		memset(TraceBack[i], 0, uPrefixCountB*sizeof(int));

// Special case for i=0
	unsigned **ptrSortOrderA = SortOrderA;
	FCOUNT **ptrFreqsA = FreqsA;
	assert(ptrSortOrderA == &(SortOrderA[0]));
	assert(ptrFreqsA == &(FreqsA[0]));
	TraceBack[0][0] = 0;

	SCORE scoreSum = 0;
	unsigned *ptrSortOrderAi = SortOrderA[0];
	const unsigned *ptrSortOrderAEnd = ptrSortOrderAi + 20;
	FCOUNT *ptrFreqsAi = FreqsA[0];
	for (; ptrSortOrderAi != ptrSortOrderAEnd; ++ptrSortOrderAi)
		{
		const unsigned uLetter = *ptrSortOrderAi;
		const FCOUNT fcLetter = ptrFreqsAi[uLetter];
		if (0 == fcLetter)
			break;
		scoreSum += fcLetter*ScoreMxB[uLetter][0];
		}
	if (0 == scoreSum)
		MPrev[0] = -2.5;
	else
		{
#if	OCC
		MPrev[0] = (logf(scoreSum) - g_scoreCenter)*OccA[0]*OccB[0];
#else
		MPrev[0] = (logf(scoreSum) - g_scoreCenter);
#endif
		}

// D(0,0) is -infinity (requires I->D).
	DPrev[0] = MINUS_INFINITY;

	for (unsigned j = 1; j < uLengthB; ++j)
		{
	// Only way to get M(0, j) looks like this:
	//		A	----X
	//		B	XXXXX
	//			0   j
	// So gap-open at j=0, gap-close at j-1.
		SCORE scoreSum = 0;
		unsigned *ptrSortOrderAi = SortOrderA[0];
		const unsigned *ptrSortOrderAEnd = ptrSortOrderAi + 20;
		FCOUNT *ptrFreqsAi = FreqsA[0];
		for (; ptrSortOrderAi != ptrSortOrderAEnd; ++ptrSortOrderAi)
			{
			const unsigned uLetter = *ptrSortOrderAi;
			const FCOUNT fcLetter = ptrFreqsAi[uLetter];
			if (0 == fcLetter)
				break;
			scoreSum += fcLetter*ScoreMxB[uLetter][j];
			}
		if (0 == scoreSum)
			MPrev[j] = -2.5;
		else
			{
#if	OCC
			MPrev[j] = (logf(scoreSum) - g_scoreCenter)*OccA[0]*OccB[j] +
			  GapOpenB[0] + GapCloseB[j-1];
#else
			MPrev[j] = (logf(scoreSum) - g_scoreCenter) +
			  GapOpenB[0] + GapCloseB[j-1];
#endif
			}
		TraceBack[0][j] = -(int) j;

	// Assume no D->I transitions, then can't be a delete if only
	// one letter from A.
		DPrev[j] = MINUS_INFINITY;
		}

	SCORE IPrev_j_1;
	for (unsigned i = 1; i < uLengthA; ++i)
		{
		++ptrSortOrderA;
		++ptrFreqsA;
		assert(ptrSortOrderA == &(SortOrderA[i]));
		assert(ptrFreqsA == &(FreqsA[i]));

		SCORE *ptrMCurr_j = MCurr;
		memset(ptrMCurr_j, 0, uLengthB*sizeof(SCORE));
		const FCOUNT *FreqsAi = *ptrFreqsA;

		const unsigned *SortOrderAi = *ptrSortOrderA;
		const unsigned *ptrSortOrderAiEnd = SortOrderAi + 20;
		const SCORE *ptrMCurrMax = MCurr + uLengthB;
		for (const unsigned *ptrSortOrderAi = SortOrderAi;
		  ptrSortOrderAi != ptrSortOrderAiEnd;
		  ++ptrSortOrderAi)
			{
			const unsigned uLetter = *ptrSortOrderAi;
			SCORE *NSBR_Letter = ScoreMxB[uLetter];
			const FCOUNT fcLetter = FreqsAi[uLetter];
			if (0 == fcLetter)
				break;
			SCORE *ptrNSBR = NSBR_Letter;
			for (SCORE *ptrMCurr = MCurr; ptrMCurr != ptrMCurrMax; ++ptrMCurr)
				*ptrMCurr += fcLetter*(*ptrNSBR++);
			}

#if	OCC
		const FCOUNT OccAi = OccA[i];
#endif
		for (unsigned j = 0; j < uLengthB; ++j)
			{
			if (MCurr[j] == 0)
				MCurr[j] = -2.5;
			else
#if	OCC
				MCurr[j] = (logf(MCurr[j]) - g_scoreCenter)*OccAi*OccB[j];
#else
				MCurr[j] = (logf(MCurr[j]) - g_scoreCenter);
#endif
			}

		ptrMCurr_j = MCurr;
		unsigned *ptrDeletePos = uDeletePos;

	// Special case for j=0
	// Only way to get M(i, 0) looks like this:
	//			0   i
	//		A	XXXXX
	//		B	----X
	// So gap-open at i=0, gap-close at i-1.
		assert(ptrMCurr_j == &(MCurr[0]));
		*ptrMCurr_j += GapOpenA[0] + GapCloseA[i-1];

		++ptrMCurr_j;

		int *ptrTraceBack_ij = TraceBack[i];
		*ptrTraceBack_ij++ = (int) i;

		SCORE *ptrMPrev_j = MPrev;
		SCORE *ptrDPrev = DPrev;
		SCORE d = *ptrDPrev;
		SCORE DNew = *ptrMPrev_j + GapOpenA[i];
		if (DNew > d)
			{
			d = DNew;
			*ptrDeletePos = i;
			}

		SCORE *ptrDCurr = DCurr;

		assert(ptrDCurr == &(DCurr[0]));
		*ptrDCurr = d;

	// Can't have an insert if no letters from B
		IPrev_j_1 = MINUS_INFINITY;

		unsigned uInsertPos = 0;
		const SCORE scoreGapOpenAi = GapOpenA[i];
		const SCORE scoreGapCloseAi_1 = GapCloseA[i-1];

		for (unsigned j = 1; j < uLengthB; ++j)
			{
		// Here, MPrev_j is preserved from previous
		// iteration so with current i,j is M[i-1][j-1]
			SCORE MPrev_j = *ptrMPrev_j;
			SCORE INew = MPrev_j + GapOpenB[j];
			if (INew > IPrev_j_1)
				{
				IPrev_j_1 = INew;
				uInsertPos = j;
				}

			SCORE scoreMax = MPrev_j;

			assert(ptrDPrev == &(DPrev[j-1]));
			SCORE scoreD = *ptrDPrev++ + scoreGapCloseAi_1;
			if (scoreD > scoreMax)
				{
				scoreMax = scoreD;
				assert(ptrDeletePos == &(uDeletePos[j-1]));
				*ptrTraceBack_ij = (int) i - (int) *ptrDeletePos;
				assert(*ptrTraceBack_ij > 0);
				}
			++ptrDeletePos;

			SCORE scoreI = IPrev_j_1 + GapCloseB[j-1];
			if (scoreI > scoreMax)
				{
				scoreMax = scoreI;
				*ptrTraceBack_ij = (int) uInsertPos - (int) j;
				assert(*ptrTraceBack_ij < 0);
				}

			assert(ptrSortOrderA == &(SortOrderA[i]));
			assert(ptrFreqsA == &(FreqsA[i]));

			*ptrMCurr_j += scoreMax;
			assert(ptrMCurr_j == &(MCurr[j]));
			++ptrMCurr_j;

			MPrev_j = *(++ptrMPrev_j);
			assert(ptrDPrev == &(DPrev[j]));
			SCORE d = *ptrDPrev;
			SCORE DNew = MPrev_j + scoreGapOpenAi;
			if (DNew > d)
				{
				d = DNew;
				assert(ptrDeletePos == &uDeletePos[j]);
				*ptrDeletePos = i;
				}
			assert(ptrDCurr + 1 == &(DCurr[j]));
			*(++ptrDCurr) = d;

			++ptrTraceBack_ij;
			}

		Rotate(MPrev, MCurr, MWork);
		Rotate(DPrev, DCurr, DWork);
		}

// Special case for i=uLengthA
	SCORE IPrev = MINUS_INFINITY;

	unsigned uInsertPos;

	for (unsigned j = 1; j < uLengthB; ++j)
		{
		SCORE INew = MPrev[j-1] + GapOpenB[j];
		if (INew > IPrev)
			{
			uInsertPos = j;
			IPrev = INew;
			}
		}

// Special case for i=uLengthA, j=uLengthB
	SCORE scoreMax = MPrev[uLengthB-1];
	int iTraceBack = 0;

	SCORE scoreD = DPrev[uLengthB-1] + GapCloseA[uLengthA-1];
	if (scoreD > scoreMax)
		{
		scoreMax = scoreD;
		iTraceBack = (int) uLengthA - (int) uDeletePos[uLengthB-1];
		}

	SCORE scoreI = IPrev + GapCloseB[uLengthB-1];
	if (scoreI > scoreMax)
		{
		scoreMax = scoreI;
		iTraceBack = (int) uInsertPos - (int) uLengthB;
		}

	TraceBack[uLengthA][uLengthB] = iTraceBack;

	TraceBackToPath(TraceBack, uLengthA, uLengthB, Path);

	return scoreMax;
	}
