#include "myutils.h"
#include "mx.h"

list<MxBase *> *MxBase::m_Matrices = 0;
unsigned MxBase::m_AllocCount;
unsigned MxBase::m_ZeroAllocCount;
unsigned MxBase::m_GrowAllocCount;
double MxBase::m_TotalBytes;
double MxBase::m_MaxBytes;

static const char *LogizeStr(const char *s)
	{
	double d = atof(s);
	d = log(d);
	return TypeToStr<float>(float(d));
	}

static const char *ExpizeStr(const char *s)
	{
	double d = atof(s);
	d = exp(d);
	return TypeToStr<float>(float(d));
	}

void MxBase::OnCtor(MxBase *Mx)
	{
#if	0
	Lock();
	if (m_Matrices == 0)
		m_Matrices = new list<MxBase *>;
	asserta(m_Matrices != 0);
	m_Matrices->push_front(Mx);
	Unlock();
#endif
	}

void MxBase::OnDtor(MxBase *Mx)
	{
#if	0
	Lock();
	if (m_Matrices == 0)
		{
		Warning("MxBase::OnDtor, m_Matrices = 0");
		Unlock();
		return;
		}
	for (list<MxBase*>::iterator p = m_Matrices->begin();
	  p != m_Matrices->end(); ++p)
		{
		if (*p == Mx)
			{
			m_Matrices->erase(p);
			if (m_Matrices->empty())
				delete m_Matrices;
			return;
			}
		}
	Warning("MxBase::OnDtor, not found");
	Unlock();
#endif // 0
	}

void MxBase::Alloc(unsigned RowCount, unsigned ColCount, const string &Name)
	{
	m_Name = Name;
	++m_AllocCount;
	if (m_AllocatedRowCount == 0)
		++m_ZeroAllocCount;

	if (RowCount > m_AllocatedRowCount || ColCount > m_AllocatedColCount)
		{
		if (m_AllocatedRowCount > 0)
			{
			++m_GrowAllocCount;
			}

		m_TotalBytes -= GetTotalBytes();

		FreeData();

		unsigned N = max(RowCount + 16, m_AllocatedRowCount);
		unsigned M = max(ColCount + 16, m_AllocatedColCount);

		AllocData(N, M);

		m_TotalBytes += GetTotalBytes();
		if (m_TotalBytes > m_MaxBytes)
			m_MaxBytes = m_TotalBytes;
		}
	
	m_Name = Name;
	asserta(RowCount <= m_AllocatedRowCount);
	asserta(ColCount <= m_AllocatedColCount);
	m_RowCount = RowCount;
	m_ColCount = ColCount;
	}

void MxBase::LogMe(bool WithData, int Opts) const
	{
	Log("\n");
	if (Opts & OPT_EXP)
		Log("Exp ");
	else if (Opts & OPT_LOG)
		Log("Log ");
	bool ZeroBased = ((Opts & OPT_ZERO_BASED) != 0);
	Log("%s(%p) Rows %u/%u, Cols %u/%u",
	  m_Name.c_str(), this,
	  m_RowCount, m_AllocatedRowCount,
	  m_ColCount, m_AllocatedColCount);
	Log("\n");
	if (!WithData || m_RowCount == 0 || m_ColCount == 0)
		return;

	const char *z = GetAsStr(0, 0);
	unsigned Width = (unsigned) strlen(z);
	unsigned Mod = 1;
	for (unsigned i = 0; i < Width; ++i)
		Mod *= 10;

	const byte *A = 0;
	const byte *B = 0;

	if (B != 0)
		{
		if (A != 0)
			Log("  ");
		Log("%5.5s", "");
		if (ZeroBased)
			for (unsigned j = 0; j < m_ColCount; ++j)
				Log("%*c", Width, B[j]);
		else
			for (unsigned j = 0; j < m_ColCount; ++j)
				Log("%*c", Width, j == 0 ? ' ' : B[j-1]);
		Log("\n");
		}

	if (A != 0)
		Log("  ");
	Log("%5.5s", "");
	for (unsigned j = 0; j < m_ColCount; ++j)
		Log("%*u", Width, j%Mod);
	Log("\n");

	for (unsigned i = 0; i < m_RowCount; ++i)
		{
		if (A != 0)
			{
			if (ZeroBased)
				Log("%c ", A[i]);
			else
				Log("%c ", i == 0 ? ' ' : A[i-1]);
			}
		Log("%4u ", i);
		
		for (unsigned j = 0; j < m_ColCount; ++j)
			{
			const char *s = GetAsStr(i, j);
			if (Opts & OPT_LOG)
				s = LogizeStr(s);
			else if (Opts & OPT_EXP)
				s = ExpizeStr(s);
			Log("%s", s);
			}
		Log("\n");
		}
	}
static unsigned g_MatrixFileCount;

void MxBase::LogCounts()
	{
	Log("\n");
	Log("MxBase::LogCounts()\n");
	Log("      What           N\n");
	Log("----------  ----------\n");
	Log("    Allocs  %10u\n", m_AllocCount);
	Log("ZeroAllocs  %10u\n", m_ZeroAllocCount);
	Log("     Grows  %10u\n", m_GrowAllocCount);
	Log("     Bytes  %10.10s\n", MemBytesToStr(m_TotalBytes));
	Log(" Max bytes  %10.10s\n", MemBytesToStr(m_MaxBytes));
	}

