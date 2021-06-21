#include "muscle.h"
#include "distfunc.h"
#include "seqvect.h"
#include <math.h>

const unsigned TRIPLE_COUNT = 20*20*20;

struct TripleCount
	{
	unsigned m_uSeqCount;			// How many sequences have this triple?
	unsigned int *m_Counts;		// m_Counts[s] = nr of times triple found in seq s
	};
static TripleCount *TripleCounts;

// WARNING: Sequences MUST be stripped of gaps and upper case!
void DistKmer20_3(const SeqVect &v, DistFunc &DF)
	{
	const unsigned uSeqCount = v.Length();

	DF.SetCount(uSeqCount);
	if (0 == uSeqCount)
		return;
	for (unsigned uSeq1 = 0; uSeq1 < uSeqCount; ++uSeq1)
		{
		DF.SetDist(uSeq1, uSeq1, 0);
		for (unsigned uSeq2 = 0; uSeq2 < uSeq1; ++uSeq2)
			DF.SetDist(uSeq1, uSeq2, 0);
		}

	const unsigned uTripleArrayBytes = TRIPLE_COUNT*sizeof(TripleCount);
	TripleCounts = (TripleCount *) malloc(uTripleArrayBytes);
	if (0 == TripleCounts)
		Quit("Not enough memory (TripleCounts)");
	memset(TripleCounts, 0, uTripleArrayBytes);

	for (unsigned uWord = 0; uWord < TRIPLE_COUNT; ++uWord)
		{
		TripleCount &tc = *(TripleCounts + uWord);
		const unsigned uBytes = uSeqCount*sizeof(int);
		tc.m_Counts = (unsigned int *) malloc(uBytes);
		memset(tc.m_Counts, 0, uBytes);
		}

	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		Seq &s = *(v[uSeqIndex]);
		const unsigned uSeqLength = s.Length();
		for (unsigned uPos = 0; uPos < uSeqLength - 2; ++uPos)
			{
			const unsigned uLetter1 = CharToLetterEx(s[uPos]);
			if (uLetter1 >= 20)
				continue;
			const unsigned uLetter2 = CharToLetterEx(s[uPos+1]);
			if (uLetter2 >= 20)
				continue;
			const unsigned uLetter3 = CharToLetterEx(s[uPos+2]);
			if (uLetter3 >= 20)
				continue;

			const unsigned uWord = uLetter1 + uLetter2*20 + uLetter3*20*20;
			assert(uWord < TRIPLE_COUNT);

			TripleCount &tc = *(TripleCounts + uWord);
			const unsigned uOldCount = tc.m_Counts[uSeqIndex];
			if (0 == uOldCount)
				++(tc.m_uSeqCount);

			++(tc.m_Counts[uSeqIndex]);
			}
		}

#if TRACE
	{
	Log("TripleCounts\n");
	unsigned uGrandTotal = 0;
	for (unsigned uWord = 0; uWord < TRIPLE_COUNT; ++uWord)
		{
		const TripleCount &tc = *(TripleCounts + uWord);
		if (0 == tc.m_uSeqCount)
			continue;

		const unsigned uLetter3 = uWord/(20*20);
		const unsigned uLetter2 = (uWord - uLetter3*20*20)/20;
		const unsigned uLetter1 = uWord%20;
		Log("Word %6u %c%c%c   %6u",
		  uWord,
		  LetterToCharAmino(uLetter1),
		  LetterToCharAmino(uLetter2),
		  LetterToCharAmino(uLetter3),
		  tc.m_uSeqCount);

		unsigned uSeqCountWithThisWord = 0;
		for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
			{
			const unsigned uCount = tc.m_Counts[uSeqIndex];
			if (uCount > 0)
				{
				++uSeqCountWithThisWord;
				Log(" %u=%u", uSeqIndex, uCount);
				uGrandTotal += uCount;
				}
			}
		if (uSeqCountWithThisWord != tc.m_uSeqCount)
			Log(" *** SQ ERROR *** %u %u", tc.m_uSeqCount, uSeqCountWithThisWord);
		Log("\n");
		}
	
	unsigned uTotalBySeqLength = 0;
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		Seq &s = *(v[uSeqIndex]);
		const unsigned uSeqLength = s.Length();
		uTotalBySeqLength += uSeqLength - 2;
		}
	if (uGrandTotal != uTotalBySeqLength)
		Log("*** TOTALS DISAGREE *** %u %u\n", uGrandTotal, uTotalBySeqLength);
	}
#endif

	const unsigned uSeqListBytes = uSeqCount*sizeof(unsigned);
	unsigned int *SeqList = (unsigned int *) malloc(uSeqListBytes);

	for (unsigned uWord = 0; uWord < TRIPLE_COUNT; ++uWord)
		{
		const TripleCount &tc = *(TripleCounts + uWord);
		if (0 == tc.m_uSeqCount)
			continue;

		unsigned uSeqCountFound = 0;
		memset(SeqList, 0, uSeqListBytes);

		for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
			{
			if (tc.m_Counts[uSeqIndex] > 0)
				{
				SeqList[uSeqCountFound] = uSeqIndex;
				++uSeqCountFound;
				if (uSeqCountFound == tc.m_uSeqCount)
					break;
				}
			}
		assert(uSeqCountFound == tc.m_uSeqCount);

		for (unsigned uSeq1 = 0; uSeq1 < uSeqCountFound; ++uSeq1)
			{
			const unsigned uSeqIndex1 = SeqList[uSeq1];
			const unsigned uCount1 = tc.m_Counts[uSeqIndex1];
			for (unsigned uSeq2 = 0; uSeq2 < uSeq1; ++uSeq2)
				{
				const unsigned uSeqIndex2 = SeqList[uSeq2];
				const unsigned uCount2 = tc.m_Counts[uSeqIndex2];
				const unsigned uMinCount = uCount1 < uCount2 ? uCount1 : uCount2;
				const double d = DF.GetDist(uSeqIndex1, uSeqIndex2);
				DF.SetDist(uSeqIndex1, uSeqIndex2, (float) (d + uMinCount));
				}
			}
		}
	delete[] SeqList;
	free(TripleCounts);

	unsigned uDone = 0;
	const unsigned uTotal = (uSeqCount*(uSeqCount - 1))/2;

	const uint PairCount = (uSeqCount*(uSeqCount - 1))/2;
	uint PairIndex = 0;
	for (unsigned uSeq1 = 0; uSeq1 < uSeqCount; ++uSeq1)
		{
		DF.SetDist(uSeq1, uSeq1, 0.0);

		const Seq &s1 = *(v[uSeq1]);
		const unsigned uLength1 = s1.Length();

		for (unsigned uSeq2 = 0; uSeq2 < uSeq1; ++uSeq2)
			{
			ProgressStep(PairIndex++, PairCount, "K-mer dist(20,3)");
			const Seq &s2 = *(v[uSeq2]);
			const unsigned uLength2 = s2.Length();
			unsigned uMinLength = uLength1 < uLength2 ? uLength1 : uLength2;
			if (uMinLength < 3)
				{
				DF.SetDist(uSeq1, uSeq2, 1.0);
				continue;
				}

			const double dTripleCount = DF.GetDist(uSeq1, uSeq2);
			if (dTripleCount == 0)
				{
				DF.SetDist(uSeq1, uSeq2, 1.0);
				continue;
				}
			double dNormalizedTripletScore = dTripleCount/(uMinLength - 2);
			DF.SetDist(uSeq1, uSeq2, (float) dNormalizedTripletScore);
#if	TRACE
			{
			Log("%s - %s  Triplet count = %g  Lengths %u, %u Estimated pwid = %g\n",
			  s1.GetName(), s2.GetName(), dTripleCount, uLength1, uLength2,
			  dEstimatedPairwiseIdentity);
			}
#endif
			}
		}
	}
