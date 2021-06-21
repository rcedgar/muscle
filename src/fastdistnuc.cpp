#include "muscle.h"
#include "distfunc.h"
#include "seqvect.h"
#include <math.h>

#define TRACE 0

#define MIN(x, y)	(((x) < (y)) ? (x) : (y))
#define MAX(x, y)	(((x) > (y)) ? (x) : (y))

const unsigned TUPLE_COUNT = 6*6*6*6*6*6;
static unsigned char Count1[TUPLE_COUNT];
static unsigned char Count2[TUPLE_COUNT];

// Nucleotide groups according to MAFFT (sextet5)
// 0 =  A
// 1 =  C
// 2 =  G
// 3 =  T
// 4 =  other

static unsigned ResidueGroup[] =
	{
	0,		// NX_A,
	1,		// NX_C,
	2,		// NX_G,
	3,		// NX_T/U
	4,		// NX_N,
	4,		// NX_R,
	4,		// NX_Y,
	4,		// NX_GAP
	};
static unsigned uResidueGroupCount = sizeof(ResidueGroup)/sizeof(ResidueGroup[0]);

static char *TupleToStr(int t)
	{
	static char s[7];
	int t1, t2, t3, t4, t5, t6;

	t1 = t%6;
	t2 = (t/6)%6;
	t3 = (t/(6*6))%6;
	t4 = (t/(6*6*6))%6;
	t5 = (t/(6*6*6*6))%6;
	t6 = (t/(6*6*6*6*6))%6;

	s[5] = '0' + t1;
	s[4] = '0' + t2;
	s[3] = '0' + t3;
	s[2] = '0' + t4;
	s[1] = '0' + t5;
	s[0] = '0' + t6;
	return s;
	}

static unsigned GetTuple(const unsigned uLetters[], unsigned n)
	{
	assert(uLetters[n] < uResidueGroupCount);
	assert(uLetters[n+1] < uResidueGroupCount);
	assert(uLetters[n+2] < uResidueGroupCount);
	assert(uLetters[n+3] < uResidueGroupCount);
	assert(uLetters[n+4] < uResidueGroupCount);
	assert(uLetters[n+5] < uResidueGroupCount);

	unsigned u1 = ResidueGroup[uLetters[n]];
	unsigned u2 = ResidueGroup[uLetters[n+1]];
	unsigned u3 = ResidueGroup[uLetters[n+2]];
	unsigned u4 = ResidueGroup[uLetters[n+3]];
	unsigned u5 = ResidueGroup[uLetters[n+4]];
	unsigned u6 = ResidueGroup[uLetters[n+5]];

	return u6 + u5*6 + u4*6*6 + u3*6*6*6 + u2*6*6*6*6 + u1*6*6*6*6*6;
	}

static void CountTuples(const unsigned L[], unsigned uTupleCount, unsigned char Count[])
	{
	memset(Count, 0, TUPLE_COUNT*sizeof(unsigned char));
	for (unsigned n = 0; n < uTupleCount; ++n)
		{
		const unsigned uTuple = GetTuple(L, n);
		++(Count[uTuple]);
		}
	}

static void ListCount(const unsigned char Count[])
	{
	for (unsigned n = 0; n < TUPLE_COUNT; ++n)
		{
		if (0 == Count[n])
			continue;
		Log("%s  %u\n", TupleToStr(n), Count[n]);
		}
	}

void DistKmer4_6(const SeqVect &v, DistFunc &DF)
	{
	if (ALPHA_DNA != g_Alpha && ALPHA_RNA != g_Alpha)
		Quit("DistKmer4_6 requires nucleo alphabet");

	const unsigned uSeqCount = v.Length();

	DF.SetCount(uSeqCount);
	if (0 == uSeqCount)
		return;

// Initialize distance matrix to zero
	for (unsigned uSeq1 = 0; uSeq1 < uSeqCount; ++uSeq1)
		{
		DF.SetDist(uSeq1, uSeq1, 0);
		for (unsigned uSeq2 = 0; uSeq2 < uSeq1; ++uSeq2)
			DF.SetDist(uSeq1, uSeq2, 0);
		}

// Convert to letters
	unsigned **Letters = new unsigned *[uSeqCount];
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		Seq &s = *(v[uSeqIndex]);
		const unsigned uSeqLength = s.Length();
		unsigned *L = new unsigned[uSeqLength];
		Letters[uSeqIndex] = L;
		for (unsigned n = 0; n < uSeqLength; ++n)
			{
			char c = s[n];
			L[n] = CharToLetterEx(c);
			if (L[n] >= 4)
				L[n] = 4;
			}
		}

	unsigned **uCommonTupleCount = new unsigned *[uSeqCount];
	for (unsigned n = 0; n < uSeqCount; ++n)
		{
		uCommonTupleCount[n] = new unsigned[uSeqCount];
		memset(uCommonTupleCount[n], 0, uSeqCount*sizeof(unsigned));
		}

	const unsigned uPairCount = (uSeqCount*(uSeqCount + 1))/2;
	unsigned uCount = 0;
	for (unsigned uSeq1 = 0; uSeq1 < uSeqCount; ++uSeq1)
		{
		Seq &seq1 = *(v[uSeq1]);
		const unsigned uSeqLength1 = seq1.Length();
		if (uSeqLength1 < 5)
			continue;

		const unsigned uTupleCount = uSeqLength1 - 5;
		const unsigned *L = Letters[uSeq1];
		CountTuples(L, uTupleCount, Count1);
#if	TRACE
		{
		Log("Seq1=%d\n", uSeq1);
		Log("Groups:\n");
		for (unsigned n = 0; n < uSeqLength1; ++n)
			Log("%u", ResidueGroup[L[n]]);
		Log("\n");

		Log("Tuples:\n");
		ListCount(Count1);
		}
#endif

		SetProgressDesc("K-mer dist pass 1");
		for (unsigned uSeq2 = 0; uSeq2 <= uSeq1; ++uSeq2)
			{
			if (0 == uCount%500)
				Progress(uCount, uPairCount);
			++uCount;
			Seq &seq2 = *(v[uSeq2]);
			const unsigned uSeqLength2 = seq2.Length();
			if (uSeqLength2 < 5)
				{
				if (uSeq1 == uSeq2)
					DF.SetDist(uSeq1, uSeq2, 0);
				else
					DF.SetDist(uSeq1, uSeq2, 1);
				continue;
				}

		// First pass through seq 2 to count tuples
			const unsigned uTupleCount = uSeqLength2 - 5;
			const unsigned *L = Letters[uSeq2];
			CountTuples(L, uTupleCount, Count2);
#if	TRACE
			Log("Seq2=%d Counts=\n", uSeq2);
			ListCount(Count2);
#endif

		// Second pass to accumulate sum of shared tuples
		// MAFFT defines this as the sum over unique tuples
		// in seq2 of the minimum of the number of tuples found
		// in the two sequences.
			unsigned uSum = 0;
			for (unsigned n = 0; n < uTupleCount; ++n)
				{
				const unsigned uTuple = GetTuple(L, n);
				uSum += MIN(Count1[uTuple], Count2[uTuple]);

			// This is a hack to make sure each unique tuple counted only once.
				Count2[uTuple] = 0;
				}
#if	TRACE
			{
			Seq &s1 = *(v[uSeq1]);
			Seq &s2 = *(v[uSeq2]);
			const char *pName1 = s1.GetName();
			const char *pName2 = s2.GetName();
			Log("Common count %s(%d) - %s(%d) =%u\n",
			  pName1, uSeq1, pName2, uSeq2, uSum);
			}
#endif
			uCommonTupleCount[uSeq1][uSeq2] = uSum;
			uCommonTupleCount[uSeq2][uSeq1] = uSum;
			}
		}
	ProgressStepsDone();

	uCount = 0;
	SetProgressDesc("K-mer dist pass 2");
	for (unsigned uSeq1 = 0; uSeq1 < uSeqCount; ++uSeq1)
		{
		Seq &s1 = *(v[uSeq1]);
		const char *pName1 = s1.GetName();

		double dCommonTupleCount11 = uCommonTupleCount[uSeq1][uSeq1];
		if (0 == dCommonTupleCount11)
			dCommonTupleCount11 = 1;

		DF.SetDist(uSeq1, uSeq1, 0);
		for (unsigned uSeq2 = 0; uSeq2 < uSeq1; ++uSeq2)
			{
			if (0 == uCount%500)
				Progress(uCount, uPairCount);
			++uCount;

			double dCommonTupleCount22 = uCommonTupleCount[uSeq2][uSeq2];
			if (0 == dCommonTupleCount22)
				dCommonTupleCount22 = 1;

			const double dDist1 = 3.0*(dCommonTupleCount11 - uCommonTupleCount[uSeq1][uSeq2])
			  /dCommonTupleCount11;
			const double dDist2 = 3.0*(dCommonTupleCount22 - uCommonTupleCount[uSeq1][uSeq2])
			  /dCommonTupleCount22;

		// dMinDist is the value used for tree-building in MAFFT
			const double dMinDist = MIN(dDist1, dDist2);
			DF.SetDist(uSeq1, uSeq2, (float) dMinDist);

			//const double dEstimatedPctId = TupleDistToEstimatedPctId(dMinDist);
			//g_dfPwId.SetDist(uSeq1, uSeq2, dEstimatedPctId);
		// **** TODO **** why does this make score slightly worse??
			//const double dKimuraDist = KimuraDist(dEstimatedPctId);
			//DF.SetDist(uSeq1, uSeq2, dKimuraDist);
			}
		}
	ProgressStepsDone();

	for (unsigned n = 0; n < uSeqCount; ++n)
		{
		delete[] uCommonTupleCount[n];
		delete[] Letters[n];
		}
	delete[] uCommonTupleCount;
	delete[] Letters;
	}
