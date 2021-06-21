#ifndef dpreglist_h
#define dpreglist_h

#include "diaglist.h"

enum DPREGIONTYPE
	{
	DPREGIONTYPE_Unknown,
	DPREGIONTYPE_Diag,
	DPREGIONTYPE_Rect
	};

struct DPRegion
	{
	DPREGIONTYPE m_Type;
	union
		{
		Diag m_Diag;
		Rect m_Rect;
		};
	};

const unsigned MAX_DPREGIONS = 1024;

class DPRegionList
	{
public:
	DPRegionList()
		{
		m_uCount = 0;
		}
	~DPRegionList()
		{
		Free();
		}

public:
// Creation
	void Clear()
		{
		Free();
		}
	void Add(const DPRegion &r);

// Accessors
	unsigned GetCount() const
		{
		return m_uCount;
		}

	const DPRegion &Get(unsigned uIndex) const
		{
		assert(uIndex < m_uCount);
		return m_DPRegions[uIndex];
		}

	unsigned GetDPArea() const;

// Diagnostics
	void LogMe() const;

private:
	void Free()
		{
		m_uCount = 0;
		}

private:
	unsigned m_uCount;
	DPRegion m_DPRegions[MAX_DPREGIONS];
	};

void DiagListToDPRegionList(const DiagList &DL, DPRegionList &RL,
  unsigned uLengthA, unsigned uLengthB);

#endif	// dpreglist_h
