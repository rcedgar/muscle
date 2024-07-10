#include "muscle.h"
#include "msa.h"
#include "tree.h"
#include "seq.h"

// These global variables are a hack to allow the tree
// dependent iteration code to communicate the edge
// used to divide the tree. The three-way weighting
// scheme needs to know this edge in order to compute
// sequence weights.
static const Tree *g_ptrMuscleTree = 0;
unsigned g_uTreeSplitNode1 = NULL_NEIGHBOR;
unsigned g_uTreeSplitNode2 = NULL_NEIGHBOR;

void MSAFromSeqRange(const MSA &msaIn, unsigned uFromSeqIndex, unsigned uSeqCount,
  MSA &msaOut)
	{
	const unsigned uColCount = msaIn.GetColCount();
	msaOut.SetSize(uSeqCount, uColCount);

	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		const char *ptrName = msaIn.GetSeqName(uFromSeqIndex + uSeqIndex);
		msaOut.SetSeqName(uSeqIndex, ptrName);

		for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
			{
			const char c = msaIn.GetChar(uFromSeqIndex + uSeqIndex, uColIndex);
			msaOut.SetChar(uSeqIndex, uColIndex, c);
			}
		}
	}

void MSAFromColRange(const MSA &msaIn, unsigned uFromColIndex, unsigned uColCount,
  MSA &msaOut)
	{
	const unsigned uSeqCount = msaIn.GetSeqCount();
	const unsigned uInColCount = msaIn.GetColCount();

	if (uFromColIndex + uColCount - 1 > uInColCount)
		Die("MSAFromColRange, out of bounds");

	msaOut.SetSize(uSeqCount, uColCount);

	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		const char *ptrName = msaIn.GetSeqName(uSeqIndex);
		unsigned uId = msaIn.GetSeqId(uSeqIndex);
		msaOut.SetSeqName(uSeqIndex, ptrName);
		msaOut.SetSeqId(uSeqIndex, uId);

		for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
			{
			const char c = msaIn.GetChar(uSeqIndex, uFromColIndex + uColIndex);
			msaOut.SetChar(uSeqIndex, uColIndex, c);
			}
		}
	}

void DeleteGappedCols(MSA &msa)
	{
	unsigned uColIndex = 0;
	for (;;)
		{
		if (uColIndex >= msa.GetColCount())
			break;
		if (msa.IsGapColumn(uColIndex))
			msa.DeleteCol(uColIndex);
		else
			++uColIndex;
		}
	}

void MSAFromSeqSubset(const MSA &msaIn, const unsigned uSeqIndexes[], unsigned uSeqCount,
  MSA &msaOut)
	{
	const unsigned uColCount = msaIn.GetColCount();
	msaOut.SetSize(uSeqCount, uColCount);
	for (unsigned uSeqIndexOut = 0; uSeqIndexOut < uSeqCount; ++uSeqIndexOut)
		{
		unsigned uSeqIndexIn = uSeqIndexes[uSeqIndexOut];
		const char *ptrName = msaIn.GetSeqName(uSeqIndexIn);
		unsigned uId = msaIn.GetSeqId(uSeqIndexIn);
		msaOut.SetSeqName(uSeqIndexOut, ptrName);
		//msaOut.SetSeqId(uSeqIndexOut, uId);
		for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
			{
			const char c = msaIn.GetChar(uSeqIndexIn, uColIndex);
			msaOut.SetChar(uSeqIndexOut, uColIndex, c);
			}
		}
	}

void AssertMSAEqIgnoreCaseAndGaps(const MSA &msa1, const MSA &msa2)
	{
	const unsigned uSeqCount1 = msa1.GetSeqCount();
	const unsigned uSeqCount2 = msa2.GetSeqCount();
	if (uSeqCount1 != uSeqCount2)
		Die("Seq count differs");

	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount1; ++uSeqIndex)
		{
		Seq seq1;
		msa1.GetSeq(uSeqIndex, seq1);

		unsigned uId = msa1.GetSeqId(uSeqIndex);
		unsigned uSeqIndex2 = msa2.GetSeqIndex(uId);

		Seq seq2;
		msa2.GetSeq(uSeqIndex2, seq2);

		if (!seq1.EqIgnoreCaseAndGaps(seq2))
			{
			Log("Input:\n");
			seq1.LogMe();
			Log("Output:\n");
			seq2.LogMe();
			Die("Seq %s differ ", msa1.GetSeqName(uSeqIndex));
			}
		}
	}

void AssertMSAEq(const MSA &msa1, const MSA &msa2)
	{
	const unsigned uSeqCount1 = msa1.GetSeqCount();
	const unsigned uSeqCount2 = msa2.GetSeqCount();
	if (uSeqCount1 != uSeqCount2)
		Die("Seq count differs");

	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount1; ++uSeqIndex)
		{
		Seq seq1;
		msa1.GetSeq(uSeqIndex, seq1);

		unsigned uId = msa1.GetSeqId(uSeqIndex);
		unsigned uSeqIndex2 = msa2.GetSeqIndex(uId);

		Seq seq2;
		msa2.GetSeq(uSeqIndex2, seq2);

		if (!seq1.Eq(seq2))
			{
			Log("Input:\n");
			seq1.LogMe();
			Log("Output:\n");
			seq2.LogMe();
			Die("Seq %s differ ", msa1.GetSeqName(uSeqIndex));
			}
		}
	}

static unsigned g_uMuscleIdCount;

#define	LOCAL_VERBOSE	0

// Append msa2 at the end of msa1
void MSAAppend(MSA &msa1, const MSA &msa2)
	{
	const unsigned uSeqCount = msa1.GetSeqCount();

	const unsigned uColCount1 = msa1.GetColCount();
	const unsigned uColCount2 = msa2.GetColCount();
	const unsigned uColCountCat = uColCount1 + uColCount2;

	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		unsigned uId = msa1.GetSeqId(uSeqIndex);
		unsigned uSeqIndex2 = msa2.GetSeqIndex(uId);
		for (unsigned uColIndex = 0; uColIndex < uColCount2; ++uColIndex)
			{
			const char c = msa2.GetChar(uSeqIndex2, uColIndex);
			msa1.SetChar(uSeqIndex, uColCount1 + uColIndex, c);
			}
		}
	}

// "Catenate" two MSAs (by bad analogy with UNIX cat command).
// msa1 and msa2 must have same sequence names, but possibly
// in a different order.
// msaCat is the combined alignment produce by appending
// sequences in msa2 to sequences in msa1.
void MSACat(const MSA &msa1, const MSA &msa2, MSA &msaCat)
	{
	const unsigned uSeqCount = msa1.GetSeqCount();

	const unsigned uColCount1 = msa1.GetColCount();
	const unsigned uColCount2 = msa2.GetColCount();
	const unsigned uColCountCat = uColCount1 + uColCount2;

	msaCat.SetSize(uSeqCount, uColCountCat);

	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		for (unsigned uColIndex = 0; uColIndex < uColCount1; ++uColIndex)
			{
			const char c = msa1.GetChar(uSeqIndex, uColIndex);
			msaCat.SetChar(uSeqIndex, uColIndex, c);
			}

		const char *ptrSeqName = msa1.GetSeqName(uSeqIndex);
		unsigned uSeqIndex2;
		msaCat.SetSeqName(uSeqIndex, ptrSeqName);
		bool bFound = msa2.GetSeqIndex(ptrSeqName, &uSeqIndex2);
		if (bFound)
			{
			for (unsigned uColIndex = 0; uColIndex < uColCount2; ++uColIndex)
				{
				const char c = msa2.GetChar(uSeqIndex2, uColIndex);
				msaCat.SetChar(uSeqIndex, uColCount1 + uColIndex, c);
				}
			}
		else
			{
			for (unsigned uColIndex = 0; uColIndex < uColCount2; ++uColIndex)
				msaCat.SetChar(uSeqIndex, uColCount1 + uColIndex, '-');
			}
		}
	}
