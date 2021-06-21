#include "muscle.h"
#include "dpreglist.h"

unsigned DPRegionList::GetDPArea() const
	{
	unsigned uArea = 0;
	for (unsigned i = 0; i < m_uCount; ++i)
		{
		const DPRegion &r = m_DPRegions[i];
		if (DPREGIONTYPE_Rect == r.m_Type)
			uArea += r.m_Rect.m_uLengthA*r.m_Rect.m_uLengthB;
		}
	return uArea;
	}

void DPRegionList::Add(const DPRegion &r)
	{
	if (m_uCount == MAX_DPREGIONS)
		Quit("DPRegionList::Add, overflow %d", m_uCount);
	m_DPRegions[m_uCount] = r;
	++m_uCount;
	}

void DPRegionList::LogMe() const
	{
	Log("DPRegionList::LogMe, count=%u\n", m_uCount);
	Log("Region  Type  StartA  StartB    EndA    EndB\n");
	Log("------  ----  ------  ------    ----    ----\n");
	for (unsigned i = 0; i < m_uCount; ++i)
		{
		const DPRegion &r = m_DPRegions[i];
		Log("%6u  ", i);
		if (DPREGIONTYPE_Diag == r.m_Type)
			Log("Diag  %6u  %6u  %6u  %6u\n",
			  r.m_Diag.m_uStartPosA,
			  r.m_Diag.m_uStartPosB,
			  r.m_Diag.m_uStartPosA + r.m_Diag.m_uLength - 1,
			  r.m_Diag.m_uStartPosB + r.m_Diag.m_uLength - 1);
		else if (DPREGIONTYPE_Rect == r.m_Type)
			Log("Rect  %6u  %6u  %6u  %6u\n",
			  r.m_Rect.m_uStartPosA,
			  r.m_Rect.m_uStartPosB,
			  r.m_Rect.m_uStartPosA + r.m_Rect.m_uLengthA - 1,
			  r.m_Rect.m_uStartPosB + r.m_Rect.m_uLengthB - 1);
		else
			Log(" *** ERROR *** Type=%u\n", r.m_Type);
		}
	}

void DiagListToDPRegionList(const DiagList &DL, DPRegionList &RL,
  unsigned uLengthA, unsigned uLengthB)
	{
	if (g_uDiagMargin > g_uMinDiagLength/2)
		Quit("Invalid parameters, diagmargin=%d must be <= 2*diaglength=%d",
		  g_uDiagMargin, g_uMinDiagLength);

	unsigned uStartPosA = 0;
	unsigned uStartPosB = 0;
	const unsigned uDiagCount = DL.GetCount();
	DPRegion r;
	for (unsigned uDiagIndex = 0; uDiagIndex < uDiagCount; ++uDiagIndex)
		{
		const Diag &d = DL.Get(uDiagIndex);
		assert(d.m_uLength >= g_uMinDiagLength);
		const unsigned uStartVertexA = d.m_uStartPosA + g_uDiagMargin - 1;
		const unsigned uStartVertexB = d.m_uStartPosB + g_uDiagMargin - 1;
		const unsigned uEndVertexA = d.m_uStartPosA + d.m_uLength - g_uDiagMargin;
		const unsigned uEndVertexB = d.m_uStartPosB + d.m_uLength - g_uDiagMargin;

		r.m_Type = DPREGIONTYPE_Rect;
		r.m_Rect.m_uStartPosA = uStartPosA;
		r.m_Rect.m_uStartPosB = uStartPosB;

		assert(uStartVertexA + 1 >= uStartPosA);
		assert(uStartVertexB + 1 >= uStartPosB);
		r.m_Rect.m_uLengthA = uStartVertexA + 1 - uStartPosA;
		r.m_Rect.m_uLengthB = uStartVertexB + 1 - uStartPosB;
		RL.Add(r);

		if (uEndVertexA > uStartVertexA + 1)
			{
			const unsigned uDiagLengthMinusCaps = uEndVertexA - uStartVertexA - 1;

			r.m_Type = DPREGIONTYPE_Diag;
			r.m_Diag.m_uStartPosA = uStartVertexA + 1;
			r.m_Diag.m_uStartPosB = uStartVertexB + 1;
			assert(uEndVertexA - uStartVertexA == uEndVertexB - uStartVertexB);
			r.m_Diag.m_uLength = uEndVertexA - uStartVertexA - 1;
			RL.Add(r);
			}

		uStartPosA = uEndVertexA;
		uStartPosB = uEndVertexB;
		}

	assert((int) uLengthA - (int) uStartPosA >= (int) g_uDiagMargin);
	assert((int) uLengthB - (int) uStartPosB >= (int) g_uDiagMargin);

	r.m_Type = DPREGIONTYPE_Rect;
	r.m_Rect.m_uStartPosA = uStartPosA;
	r.m_Rect.m_uStartPosB = uStartPosB;

	assert(uLengthA >= uStartPosA);
	assert(uLengthB >= uStartPosB);
	r.m_Rect.m_uLengthA = uLengthA - uStartPosA;
	r.m_Rect.m_uLengthB = uLengthB - uStartPosB;
	RL.Add(r);
	}
