#include "muscle.h"
#include "profile.h"
#include "pwpath.h"
#include "seq.h"

extern SCOREMATRIX VTML_SP;

// #define SUBST(i, j)	Subst(seqA, seqB, i, j)
#define SUBST(i, j)		MxRowA[i][seqB.GetLetter(j)]

static SCORE Subst(const Seq &seqA, const Seq &seqB, unsigned i, unsigned j)
	{
	assert(i < seqA.Length());
	assert(j < seqB.Length());

	unsigned uLetterA = seqA.GetLetter(i);
	unsigned uLetterB = seqB.GetLetter(j);
	return VTML_SP[uLetterA][uLetterB] + g_scoreCenter;
	}

struct DP_MEMORY
	{
	unsigned uLength;
	SCORE *MPrev;
	SCORE *MCurr;
	SCORE *MWork;
	SCORE *DPrev;
	SCORE *DCurr;
	SCORE *DWork;
	SCORE **MxRowA;
	unsigned *LettersB;
	unsigned *uDeletePos;
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
			delete[] DPM.TraceBack[i];

		delete[] DPM.MPrev;
		delete[] DPM.MCurr;
		delete[] DPM.MWork;
		delete[] DPM.DPrev;
		delete[] DPM.DCurr;
		delete[] DPM.DWork;
		delete[] DPM.MxRowA;
		delete[] DPM.LettersB;
		delete[] DPM.uDeletePos;
		delete[] DPM.TraceBack;
		}

	DPM.uLength = uLength;

	DPM.MPrev = new SCORE[uLength];
	DPM.MCurr = new SCORE[uLength];
	DPM.MWork = new SCORE[uLength];

	DPM.DPrev = new SCORE[uLength];
	DPM.DCurr = new SCORE[uLength];
	DPM.DWork = new SCORE[uLength];
	DPM.MxRowA = new SCORE *[uLength];
	DPM.LettersB = new unsigned[uLength];
	DPM.uDeletePos = new unsigned[uLength];

	DPM.TraceBack = new int*[uLength];

	for (unsigned i = 0; i < uLength; ++i)
		DPM.TraceBack[i] = new int[uLength];
	}

static void RowFromSeq(const Seq &s, SCORE *Row[])
	{
	const unsigned uLength = s.Length();
	for (unsigned i = 0; i < uLength; ++i)
		{
		char c = s.GetChar(i);
		unsigned uLetter = CharToLetter(c);
		if (uLetter < 20)
			Row[i] = VTML_SP[uLetter];
		else
			Row[i] = VTML_SP[AX_X];
		}
	}

static void LettersFromSeq(const Seq &s, unsigned Letters[])
	{
	const unsigned uLength = s.Length();
	for (unsigned i = 0; i < uLength; ++i)
		{
		char c = s.GetChar(i);
		unsigned uLetter = CharToLetter(c);
		if (uLetter < 20)
			Letters[i] = uLetter;
		else
			Letters[i] = AX_X;
		}
	}

SCORE GlobalAlignSS(const Seq &seqA, const Seq &seqB, PWPath &Path)
	{
	const unsigned uLengthA = seqA.Length();
	const unsigned uLengthB = seqB.Length();
	const unsigned uPrefixCountA = uLengthA + 1;
	const unsigned uPrefixCountB = uLengthB + 1;

	AllocDPMem(uLengthA, uLengthB);

	SCORE *MPrev = DPM.MPrev;
	SCORE *MCurr = DPM.MCurr;
	SCORE *MWork = DPM.MWork;

	SCORE *DPrev = DPM.DPrev;
	SCORE *DCurr = DPM.DCurr;
	SCORE *DWork = DPM.DWork;
	SCORE **MxRowA = DPM.MxRowA;
	unsigned *LettersB = DPM.LettersB;

	RowFromSeq(seqA, MxRowA);
	LettersFromSeq(seqB, LettersB);

	unsigned *uDeletePos = DPM.uDeletePos;

	int **TraceBack = DPM.TraceBack;

#if	DEBUG
	for (unsigned i = 0; i < uPrefixCountA; ++i)
		memset(TraceBack[i], 0, uPrefixCountB*sizeof(int));
#endif

// Special case for i=0
	TraceBack[0][0] = 0;
	MPrev[0] = MxRowA[0][LettersB[0]];

// D(0,0) is -infinity (requires I->D).
	DPrev[0] = MINUS_INFINITY;

	for (unsigned j = 1; j < uLengthB; ++j)
		{
		unsigned uLetterB = LettersB[j];

	// Only way to get M(0, j) looks like this:
	//		A	----X
	//		B	XXXXX
	//			0   j
	// So gap-open at j=0, gap-close at j-1.
		MPrev[j] = MxRowA[0][uLetterB] + g_scoreGapOpen/2; // term gaps half
		TraceBack[0][j] = -(int) j;

	// Assume no D->I transitions, then can't be a delete if only
	// one letter from A.
		DPrev[j] = MINUS_INFINITY;
		}

	SCORE IPrev_j_1;
	for (unsigned i = 1; i < uLengthA; ++i)
		{
		SCORE *ptrMCurr_j = MCurr;
		memset(ptrMCurr_j, 0, uLengthB*sizeof(SCORE));

		const SCORE *RowA = MxRowA[i];
		const SCORE *ptrRowA = MxRowA[i];
		const SCORE *ptrMCurrEnd = ptrMCurr_j + uLengthB;
		unsigned *ptrLettersB = LettersB;
		for (; ptrMCurr_j != ptrMCurrEnd; ++ptrMCurr_j)
			{
			*ptrMCurr_j = RowA[*ptrLettersB];
			++ptrLettersB;
			}

		unsigned *ptrDeletePos = uDeletePos;

	// Special case for j=0
	// Only way to get M(i, 0) looks like this:
	//			0   i
	//		A	XXXXX
	//		B	----X
	// So gap-open at i=0, gap-close at i-1.
		ptrMCurr_j = MCurr;
		assert(ptrMCurr_j == &(MCurr[0]));
		*ptrMCurr_j += g_scoreGapOpen/2;	// term gaps half

		++ptrMCurr_j;

		int *ptrTraceBack_ij = TraceBack[i];
		*ptrTraceBack_ij++ = (int) i;

		SCORE *ptrMPrev_j = MPrev;
		SCORE *ptrDPrev = DPrev;
		SCORE d = *ptrDPrev;
		SCORE DNew = *ptrMPrev_j + g_scoreGapOpen;
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

		unsigned uInsertPos;

		for (unsigned j = 1; j < uLengthB; ++j)
			{
		// Here, MPrev_j is preserved from previous
		// iteration so with current i,j is M[i-1][j-1]
			SCORE MPrev_j = *ptrMPrev_j;
			SCORE INew = MPrev_j + g_scoreGapOpen;
			if (INew > IPrev_j_1)
				{
				IPrev_j_1 = INew;
				uInsertPos = j;
				}

			SCORE scoreMax = MPrev_j;

			assert(ptrDPrev == &(DPrev[j-1]));
			SCORE scoreD = *ptrDPrev++;
			if (scoreD > scoreMax)
				{
				scoreMax = scoreD;
				assert(ptrDeletePos == &(uDeletePos[j-1]));
				*ptrTraceBack_ij = (int) i - (int) *ptrDeletePos;
				assert(*ptrTraceBack_ij > 0);
				}
			++ptrDeletePos;

			SCORE scoreI = IPrev_j_1;
			if (scoreI > scoreMax)
				{
				scoreMax = scoreI;
				*ptrTraceBack_ij = (int) uInsertPos - (int) j;
				assert(*ptrTraceBack_ij < 0);
				}

			*ptrMCurr_j += scoreMax;
			assert(ptrMCurr_j == &(MCurr[j]));
			++ptrMCurr_j;

			MPrev_j = *(++ptrMPrev_j);
			assert(ptrDPrev == &(DPrev[j]));
			SCORE d = *ptrDPrev;
			SCORE DNew = MPrev_j + g_scoreGapOpen;
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
		SCORE INew = MPrev[j-1];
		if (INew > IPrev)
			{
			uInsertPos = j;
			IPrev = INew;
			}
		}

// Special case for i=uLengthA, j=uLengthB
	SCORE scoreMax = MPrev[uLengthB-1];
	int iTraceBack = 0;

	SCORE scoreD = DPrev[uLengthB-1] - g_scoreGapOpen/2;	// term gaps half
	if (scoreD > scoreMax)
		{
		scoreMax = scoreD;
		iTraceBack = (int) uLengthA - (int) uDeletePos[uLengthB-1];
		}

	SCORE scoreI = IPrev - g_scoreGapOpen/2;
	if (scoreI > scoreMax)
		{
		scoreMax = scoreI;
		iTraceBack = (int) uInsertPos - (int) uLengthB;
		}

	TraceBack[uLengthA][uLengthB] = iTraceBack;

	TraceBackToPath(TraceBack, uLengthA, uLengthB, Path);

	return scoreMax;
	}
