#include "muscle.h"
#include "diaglist.h"
#include "pwpath.h"

#define MAX(x, y)	((x) > (y) ? (x) : (y))
#define MIN(x, y)	((x) < (y) ? (x) : (y))

void DiagList::Add(const Diag &d)
	{
	if (m_uCount == MAX_DIAGS)
		Quit("DiagList::Add, overflow %u", m_uCount);
	m_Diags[m_uCount] = d;
	++m_uCount;
	}

void DiagList::Add(unsigned uStartPosA, unsigned uStartPosB, unsigned uLength)
	{
	Diag d;
	d.m_uStartPosA = uStartPosA;
	d.m_uStartPosB = uStartPosB;
	d.m_uLength = uLength;
	Add(d);
	}

const Diag &DiagList::Get(unsigned uIndex) const
	{
	if (uIndex >= m_uCount)
		Quit("DiagList::Get(%u), count=%u", uIndex, m_uCount);
	return m_Diags[uIndex];
	}

void DiagList::LogMe() const
	{
	Log("DiagList::LogMe, count=%u\n", m_uCount);
	Log("  n  StartA  StartB  Length\n");
	Log("---  ------  ------  ------\n");
	for (unsigned n = 0; n < m_uCount; ++n)
		{
		const Diag &d = m_Diags[n];
		Log("%3u  %6u  %6u  %6u\n",
		  n, d.m_uStartPosA, d.m_uStartPosB, d.m_uLength);
		}
	}

void DiagList::FromPath(const PWPath &Path)
	{
	Clear();

	const unsigned uEdgeCount = Path.GetEdgeCount();
	unsigned uLength = 0;
	unsigned uStartPosA;
	unsigned uStartPosB;
	for (unsigned uEdgeIndex = 0; uEdgeIndex < uEdgeCount; ++uEdgeIndex)
		{
		const PWEdge &Edge = Path.GetEdge(uEdgeIndex);

	// Typical cases
		if (Edge.cType == 'M')
			{
			if (0 == uLength)
				{
				uStartPosA = Edge.uPrefixLengthA - 1;
				uStartPosB = Edge.uPrefixLengthB - 1;
				}
			++uLength;
			}
		else
			{
			if (uLength >= g_uMinDiagLength)
				Add(uStartPosA, uStartPosB, uLength);
			uLength = 0;
			}
		}

// Special case for last edge
	if (uLength >= g_uMinDiagLength)
		Add(uStartPosA, uStartPosB, uLength);
	}

bool DiagList::NonZeroIntersection(const Diag &d) const
	{
	for (unsigned n = 0; n < m_uCount; ++n)
		{
		const Diag &d2 = m_Diags[n];
		if (DiagOverlap(d, d2) > 0)
			return true;
		}
	return false;
	}

// DialogOverlap returns the length of the overlapping
// section of the two diagonals along the diagonals
// themselves; in other words, the length of
// the intersection of the two sets of cells in
// the matrix.
unsigned DiagOverlap(const Diag &d1, const Diag &d2)
	{
// Determine where the diagonals intersect the A
// axis (extending them if required). If they
// intersect at different points, they do not
// overlap. Coordinates on a diagonal are
// given by B = A + c where c is the value of
// A at the intersection with the A axis.
// Hence, c = B - A for any point on the diagonal.
	int c1 = (int) d1.m_uStartPosB - (int) d1.m_uStartPosA;
	int c2 = (int) d2.m_uStartPosB - (int) d2.m_uStartPosA;
	if (c1 != c2)
		return 0;

	assert(DiagOverlapA(d1, d2) == DiagOverlapB(d1, d2));
	return DiagOverlapA(d1, d2);
	}

// DialogOverlapA returns the length of the overlapping
// section of the projection of the two diagonals onto
// the A axis.
unsigned DiagOverlapA(const Diag &d1, const Diag &d2)
	{
	unsigned uMaxStart = MAX(d1.m_uStartPosA, d2.m_uStartPosA);
	unsigned uMinEnd = MIN(d1.m_uStartPosA + d1.m_uLength - 1,
	  d2.m_uStartPosA + d2.m_uLength - 1);

	int iLength = (int) uMinEnd - (int) uMaxStart + 1;
	if (iLength < 0)
		return 0;
	return (unsigned) iLength;
	}

// DialogOverlapB returns the length of the overlapping
// section of the projection of the two diagonals onto
// the B axis.
unsigned DiagOverlapB(const Diag &d1, const Diag &d2)
	{
	unsigned uMaxStart = MAX(d1.m_uStartPosB, d2.m_uStartPosB);
	unsigned uMinEnd = MIN(d1.m_uStartPosB + d1.m_uLength - 1,
	  d2.m_uStartPosB + d2.m_uLength - 1);

	int iLength = (int) uMinEnd - (int) uMaxStart + 1;
	if (iLength < 0)
		return 0;
	return (unsigned) iLength;
	}

// Returns true if the two diagonals can be on the
// same path through the DP matrix. If DiagCompatible
// returns false, they cannot be in the same path
// and hence "contradict" each other.
bool DiagCompatible(const Diag &d1, const Diag &d2)
	{
	if (DiagOverlap(d1, d2) > 0)
		return true;
	return 0 == DiagOverlapA(d1, d2) && 0 == DiagOverlapB(d1, d2);
	}

// Returns the length of the "break" between two diagonals.
unsigned DiagBreak(const Diag &d1, const Diag &d2)
	{
	int c1 = (int) d1.m_uStartPosB - (int) d1.m_uStartPosA;
	int c2 = (int) d2.m_uStartPosB - (int) d2.m_uStartPosA;
	if (c1 != c2)
		return 0;

	int iMaxStart = MAX(d1.m_uStartPosA, d2.m_uStartPosA);
	int iMinEnd = MIN(d1.m_uStartPosA + d1.m_uLength - 1,
	  d2.m_uStartPosA + d1.m_uLength - 1);
	int iBreak = iMaxStart - iMinEnd - 1;
	if (iBreak < 0)
		return 0;
	return (unsigned) iBreak;
	}

// Merge diagonals that are continuations of each other with
// int breaks of up to length g_uMaxDiagBreak.
// In a sorted list of diagonals, we only have to check
// consecutive entries.
void MergeDiags(DiagList &DL)
	{
	return;
#if	DEBUG
	if (!DL.IsSorted())
		Quit("MergeDiags: !IsSorted");
#endif

// TODO: Fix this!
// Breaks must be with no offset (no gaps)
	const unsigned uCount = DL.GetCount();
	if (uCount <= 1)
		return;

	DiagList NewList;

	Diag MergedDiag;
	const Diag *ptrPrev = &DL.Get(0);
	for (unsigned i = 1; i < uCount; ++i)
		{
		const Diag *ptrDiag = &DL.Get(i);
		unsigned uBreakLength = DiagBreak(*ptrPrev, *ptrDiag);
		if (uBreakLength <= g_uMaxDiagBreak)
			{
			MergedDiag.m_uStartPosA = ptrPrev->m_uStartPosA;
			MergedDiag.m_uStartPosB = ptrPrev->m_uStartPosB;
			MergedDiag.m_uLength = ptrPrev->m_uLength + ptrDiag->m_uLength
			  + uBreakLength;
			ptrPrev = &MergedDiag;
			}
		else
			{
			NewList.Add(*ptrPrev);
			ptrPrev = ptrDiag;
			}
		}
	NewList.Add(*ptrPrev);
	DL.Copy(NewList);
	}

void DiagList::DeleteIncompatible()
	{
	assert(IsSorted());

	if (m_uCount < 2)
		return;

	bool *bFlagForDeletion = new bool[m_uCount];
	for (unsigned i = 0; i < m_uCount; ++i)
		bFlagForDeletion[i] = false;

	for (unsigned i = 0; i < m_uCount; ++i)
		{
		const Diag &di = m_Diags[i];
		for (unsigned j = i + 1; j < m_uCount; ++j)
			{
			const Diag &dj = m_Diags[j];

		// Verify sorted correctly
			assert(di.m_uStartPosA <= dj.m_uStartPosA);

		// If two diagonals are incompatible and
		// one is is much longer than the other,
		// keep the longer one.
			if (!DiagCompatible(di, dj))
				{
				if (di.m_uLength > dj.m_uLength*4)
					bFlagForDeletion[j] = true;
				else if (dj.m_uLength > di.m_uLength*4)
					bFlagForDeletion[i] = true;
				else
					{
					bFlagForDeletion[i] = true;
					bFlagForDeletion[j] = true;
					}
				}
			}
		}

	for (unsigned i = 0; i < m_uCount; ++i)
		{
		const Diag &di = m_Diags[i];
		if (bFlagForDeletion[i])
			continue;

		for (unsigned j = i + 1; j < m_uCount; ++j)
			{
			const Diag &dj = m_Diags[j];
			if (bFlagForDeletion[j])
				continue;

		// Verify sorted correctly
			assert(di.m_uStartPosA <= dj.m_uStartPosA);

		// If sort order in B different from sorted order in A,
		// either diags are incompatible or we detected a repeat
		// or permutation.
			if (di.m_uStartPosB >= dj.m_uStartPosB || !DiagCompatible(di, dj))
				{
				bFlagForDeletion[i] = true;
				bFlagForDeletion[j] = true;
				}
			}
		}

	unsigned uNewCount = 0;
	Diag *NewDiags = new Diag[m_uCount];
	for (unsigned i = 0; i < m_uCount; ++i)
		{
		if (bFlagForDeletion[i])
			continue;

		const Diag &d = m_Diags[i];
		NewDiags[uNewCount] = d;
		++uNewCount;
		}
	memcpy(m_Diags, NewDiags, uNewCount*sizeof(Diag));
	m_uCount = uNewCount;
	delete[] NewDiags;
	}

void DiagList::Copy(const DiagList &DL)
	{
	Clear();
	unsigned uCount = DL.GetCount();
	for (unsigned i = 0; i < uCount; ++i)
		Add(DL.Get(i));
	}

// Check if sorted in increasing order of m_uStartPosA
bool DiagList::IsSorted() const
	{
	return true;
	unsigned uCount = GetCount();
	for (unsigned i = 1; i < uCount; ++i)
		if (m_Diags[i-1].m_uStartPosA > m_Diags[i].m_uStartPosA)
			return false;
	return true;
	}

// Sort in increasing order of m_uStartPosA
// Dumb bubble sort, but don't care about speed
// because don't get long lists.
void DiagList::Sort()
	{
	if (m_uCount < 2)
		return;

	bool bContinue = true;
	while (bContinue)
		{
		bContinue = false;
		for (unsigned i = 0; i < m_uCount - 1; ++i)
			{
			if (m_Diags[i].m_uStartPosA > m_Diags[i+1].m_uStartPosA)
				{
				Diag Tmp = m_Diags[i];
				m_Diags[i] = m_Diags[i+1];
				m_Diags[i+1] = Tmp;
				bContinue = true;
				}
			}
		}
	}

//void TestDiag()
//	{
//	Diag d1;
//	Diag d2;
//	Diag d3;
//
//	d1.m_uStartPosA = 0;
//	d1.m_uStartPosB = 1;
//	d1.m_uLength = 32;
//
//	d2.m_uStartPosA = 55;
//	d2.m_uStartPosB = 70;
//	d2.m_uLength = 36;
//
//	d3.m_uStartPosA = 102;
//	d3.m_uStartPosB = 122;
//	d3.m_uLength = 50;
//
//	DiagList DL;
//	DL.Add(d1);
//	DL.Add(d2);
//	DL.Add(d3);
//
//	Log("Before DeleteIncompatible:\n");
//	DL.LogMe();
//	DL.DeleteIncompatible();
//
//	Log("After DeleteIncompatible:\n");
//	DL.LogMe();
//
//	MergeDiags(DL);
//	Log("After Merge:\n");
//	DL.LogMe();
//
//	DPRegionList RL;
//	DiagListToDPRegionList(DL, RL, 200, 200);
//	RL.LogMe();
//	}
