#include "myutils.h"
#include "pathinfo.h"

#define TRACE	0

void PathInfo::FreeIfBig()
	{
	if (m_BufferBytes < 4*g_MaxL)
		return;

#if	TRACE
	Log("PathInfo::FreeIfBig(this=0x%lx) m_BufferBytes=%u\n",
	  (long) this, m_BufferBytes);
#endif
	myfree(m_Buffer);
	m_Buffer = 0;
	m_BufferBytes = 0;
	}

void PathInfo::Alloc(unsigned Bytes)
	{
	FreeIfBig();
	if (Bytes < m_BufferBytes)
		return;

#if	TRACE
	Log("PathInfo::Alloc(this=0x%lx, Bytes=%u) m_BufferBytes=%u\n",
	  (long) this, Bytes, m_BufferBytes);
#endif
	myfree(m_Buffer);
	m_BufferBytes = Bytes + 128;
	m_Buffer = myalloc(char, m_BufferBytes);
	}

void PathInfo::Realloc(unsigned Bytes)
	{
	asserta(Bytes > m_BufferBytes);
	char *NewBuffer = myalloc(char, Bytes);
	if (m_ColCount > 0)
		{
		asserta(m_ColCount < m_BufferBytes);
		memcpy(NewBuffer, m_Buffer, m_ColCount + 1);
		}
	myfree(m_Buffer);
	m_Buffer = NewBuffer;
	m_BufferBytes = Bytes;
	}

void PathInfo::Alloc2(unsigned LA, unsigned LB)
	{
	FreeIfBig();
	unsigned Bytes = LA + LB + 1;
	if (Bytes > m_BufferBytes)
		Alloc(Bytes);
	}

void PathInfo::SetEmpty()
	{
	FreeIfBig();

// Prevent crash if never alloc'd
	if (m_Buffer == 0)
		Alloc(g_MaxL);
	m_Buffer[0] = 0;
	m_ColCount = 0;
	}

unsigned PathInfo::GetCounts(unsigned &M, unsigned &D, unsigned &I) const
	{
	M = 0;
	D = 0;
	I = 0;
	const char *Path = GetPath();
	for (const char *p = Path; *p; ++p)
		{
		char c = *p;
		if (c == 'M')
			++M;
		else if (c == 'D')
			++D;
		else if (c == 'I')
			++I;
		else
			asserta(false);
		}
	return M + D + I;
	}

void PathInfo::Reverse()
	{
	unsigned L = GetColCount();
	for (unsigned i = 0; i < L/2; ++i)
		swap(m_Buffer[i], m_Buffer[L-i-1]);
	}

void PathInfo::AppendPath(const PathInfo &PI)
	{
	unsigned n = PI.GetColCount();
	if (n == 0)
		return;
	Grow(n);
	const char *Path = PI.GetPath();
	memcpy(m_Buffer + m_ColCount, Path, n);
	m_ColCount += n;
	m_Buffer[m_ColCount] = 0;
	}

void PathInfo::PrependPath(const PathInfo &PI)
	{
	unsigned n = PI.GetColCount();
	Grow(n);

	unsigned ColCount = PI.GetColCount();
	if (m_ColCount > 0)
		memmove(m_Buffer + ColCount, m_Buffer, m_ColCount);

	if (ColCount > 0)
		{
		const char *Path = PI.GetPath();
		memcpy(m_Buffer, Path, ColCount);
		m_ColCount += ColCount;
		asserta(m_ColCount < m_BufferBytes);
		m_Buffer[m_ColCount] = 0;
		assert(strlen(m_Buffer) == m_ColCount);
		}
	}

void PathInfo::AppendChar(char c)
	{
	Grow(1);
	m_Buffer[m_ColCount++] = c;
	m_Buffer[m_ColCount] = 0;
	}

void PathInfo::AppendMs(unsigned Count)
	{
	Grow(Count);
	for (unsigned i = 0; i < Count; ++i)
		m_Buffer[m_ColCount++] = 'M';
	m_Buffer[m_ColCount] = 0;
	}

void PathInfo::AppendDs(unsigned Count)
	{
	for (unsigned i = 0; i < Count; ++i)
		m_Buffer[m_ColCount++] = 'D';
	m_Buffer[m_ColCount] = 0;
	}

void PathInfo::AppendIs(unsigned Count)
	{
	for (unsigned i = 0; i < Count; ++i)
		m_Buffer[m_ColCount++] = 'I';
	m_Buffer[m_ColCount] = 0;
	}

unsigned PathInfo::GetLeftICount()
	{
	for (unsigned i = 0; i < m_ColCount; ++i)
		if (m_Buffer[i] != 'I')
			return i;
	asserta(false);
	return UINT_MAX;
	}

unsigned PathInfo::TrimLeftIs()
	{
	unsigned LeftICount = GetLeftICount();
	m_ColCount -= LeftICount;
	memmove(m_Buffer, m_Buffer+LeftICount, m_ColCount);
	m_Buffer[m_ColCount] = 0;
	return LeftICount;
	}

void PathInfo::TrimRightIs()
	{
	for (unsigned i = m_ColCount - 1; i > 0; --i)
		{
		if (m_Buffer[i] == 'I')
			{
			m_Buffer[i] = 0;
			--m_ColCount;
			}
		else
			{
			assert(strlen(m_Buffer) == m_ColCount);
			return;
			}
		}
	}

void PathInfo::ToOps(string &Ops, vector<uint> &Lengths) const
	{
	Ops.clear();
	Lengths.clear();
	const char *p = GetPath();
	char LastOp = *p++;
	uint n = 1;
	for (; *p; ++p)
		{
		char Op = *p;
		if (Op == LastOp)
			++n;
		else
			{
			assert(n > 0);
			Ops.push_back(LastOp);
			Lengths.push_back(n);
			n = 1;
			LastOp = Op;
			}
		}
	assert(n > 0);
	Ops.push_back(LastOp);
	Lengths.push_back(n);

	const uint M = SIZE(Ops);
	for (uint i = 0; i < M; ++i)
		{
		char c = Ops[i];
		if (c == 'D')
			Ops[i] = 'I';
		else if (c == 'I')
			Ops[i] = 'D';
		}
	}
