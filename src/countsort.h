#ifndef countsort_h
#define countsort_h

#include "gobuff.h"

class CountSortMem
	{
public:
	static const unsigned NVEC = 8;

public:
	unsigned *m_Vecs[NVEC];
	unsigned m_VecPos[NVEC];
	unsigned m_MaxValueCount;

	GoBuff<unsigned> m_Sizes;
	GoBuff<unsigned> m_Offsets;

public:
	CountSortMem()
		{
		m_MaxValueCount = 0;
		memset_zero(m_Vecs, NVEC);
		}

	void Free()
		{
		for (unsigned i = 0; i < NVEC; ++i)
			{
			myfree(m_Vecs[i]);
			m_Vecs[i] = 0;
			}
		m_MaxValueCount = 0;
		}

	void Alloc(unsigned ValueCount)
		{
		if (ValueCount <= m_MaxValueCount)
			return;

		Free();
		
		m_MaxValueCount = ValueCount;
		for (unsigned i = 0; i < NVEC; ++i)
			m_Vecs[i] = myalloc(unsigned, m_MaxValueCount);
		}
	};

unsigned CountSortOrderDesc(const unsigned *Values, unsigned ValueCount,
  CountSortMem &Mem, unsigned *Order);
unsigned CountSortSubsetDesc(const unsigned *Values, unsigned ValueCount,
  CountSortMem &Mem, const unsigned *Subset, unsigned *Result);

#endif // countsort_h
