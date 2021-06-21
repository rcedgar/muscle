#include "muscle.h"
#include "profile.h"
#include "msa.h"
#include "pwpath.h"
#include "seqvect.h"
#include "clust.h"
#include "tree.h"

#define TRACE 0

struct Range
	{
	unsigned m_uBestColLeft;
	unsigned m_uBestColRight;
	};

static void ListVertSavings(unsigned uColCount, unsigned uAnchorColCount,
  const Range *Ranges, unsigned uRangeCount)
	{
	if (!g_bVerbose || !g_bAnchors)
		return;
	double dTotalArea = uColCount*uColCount;
	double dArea = 0.0;
	for (unsigned i = 0; i < uRangeCount; ++i)
		{
		unsigned uLength = Ranges[i].m_uBestColRight - Ranges[i].m_uBestColLeft;
		dArea += uLength*uLength;
		}
	double dPct = (dTotalArea - dArea)*100.0/dTotalArea;
	Log("Anchor columns found       %u\n", uAnchorColCount);
	Log("DP area saved by anchors   %-4.1f%%\n", dPct);
	}

static void ColsToRanges(const unsigned BestCols[], unsigned uBestColCount,
  unsigned uColCount, Range Ranges[])
	{
// N best columns produces N+1 vertical blocks.
	const unsigned uRangeCount = uBestColCount + 1;
	for (unsigned uIndex = 0; uIndex < uRangeCount ; ++uIndex)
		{
		unsigned uBestColLeft = 0;
		if (uIndex > 0)
			uBestColLeft = BestCols[uIndex-1];
		
		unsigned uBestColRight = uColCount;
		if (uIndex < uBestColCount)
			uBestColRight = BestCols[uIndex];

		Ranges[uIndex].m_uBestColLeft = uBestColLeft;
		Ranges[uIndex].m_uBestColRight = uBestColRight;
		}
	}

// Return true if any changes made
bool RefineVert(MSA &msaIn, const Tree &tree, unsigned uIters)
	{
	bool bAnyChanges = false;

	const unsigned uColCountIn = msaIn.GetColCount();
	const unsigned uSeqCountIn = msaIn.GetSeqCount();

	if (uColCountIn < 3 || uSeqCountIn < 3)
		return false;

	unsigned *AnchorCols = new unsigned[uColCountIn];
	unsigned uAnchorColCount;
	SetMSAWeightsMuscle(msaIn);
	FindAnchorCols(msaIn, AnchorCols, &uAnchorColCount);

	const unsigned uRangeCount = uAnchorColCount + 1;
	Range *Ranges = new Range[uRangeCount];

#if	TRACE
	Log("%u ranges\n", uRangeCount);
#endif

	ColsToRanges(AnchorCols, uAnchorColCount, uColCountIn, Ranges);
	ListVertSavings(uColCountIn, uAnchorColCount, Ranges, uRangeCount);

#if	TRACE
	{
	Log("Anchor cols: ");
	for (unsigned i = 0; i < uAnchorColCount; ++i)
		Log(" %u", AnchorCols[i]);
	Log("\n");

	Log("Ranges:\n");
	for (unsigned i = 0; i < uRangeCount; ++i)
		Log("%4u - %4u\n", Ranges[i].m_uBestColLeft, Ranges[i].m_uBestColRight);
	}
#endif

	delete[] AnchorCols;

	MSA msaOut;
	msaOut.SetSize(uSeqCountIn, 0);

	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCountIn; ++uSeqIndex)
		{
		const char *ptrName = msaIn.GetSeqName(uSeqIndex);
		unsigned uId = msaIn.GetSeqId(uSeqIndex);
		msaOut.SetSeqName(uSeqIndex, ptrName);
		msaOut.SetSeqId(uSeqIndex, uId);
		}

	for (unsigned uRangeIndex = 0; uRangeIndex < uRangeCount; ++uRangeIndex)
		{
		MSA msaRange;

		const Range &r = Ranges[uRangeIndex];

		const unsigned uFromColIndex = r.m_uBestColLeft;
		const unsigned uRangeColCount = r.m_uBestColRight - uFromColIndex;

		if (0 == uRangeColCount)
			continue;
		else if (1 == uRangeColCount)
			{
			MSAFromColRange(msaIn, uFromColIndex, 1, msaRange);
			MSAAppend(msaOut, msaRange);
			continue;
			}
		MSAFromColRange(msaIn, uFromColIndex, uRangeColCount, msaRange);

#if	TRACE
		Log("\n-------------\n");
		Log("Range %u - %u count=%u\n", r.m_uBestColLeft, r.m_uBestColRight, uRangeColCount);
		Log("Before:\n");
		msaRange.LogMe();
#endif

		bool bLockLeft = (0 != uRangeIndex);
		bool bLockRight = (uRangeCount - 1 != uRangeIndex);
		bool bAnyChangesThisBlock = RefineHoriz(msaRange, tree, uIters, bLockLeft, bLockRight);
		bAnyChanges = (bAnyChanges || bAnyChangesThisBlock);

#if	TRACE
		Log("After:\n");
		msaRange.LogMe();
#endif

		MSAAppend(msaOut, msaRange);

#if	TRACE
		Log("msaOut after Cat:\n");
		msaOut.LogMe();
#endif
		}

#if	DEBUG
// Sanity check
	AssertMSAEqIgnoreCaseAndGaps(msaIn, msaOut);
#endif

	delete[] Ranges;
	if (bAnyChanges)
		msaIn.Copy(msaOut);
	return bAnyChanges;
	}
