#ifndef xdpmem_h
#define xdpmem_h

#include "mx.h"

static const unsigned g_MaxL = 4096;
static const double g_MaxL2 = double(g_MaxL)*double(g_MaxL);

unsigned GetSubL(unsigned L);

class ObjMgr;

class XDPMem
	{
public:
	unsigned m_MaxLA;
	unsigned m_MaxLB;
	Mx<byte> m_TBBit;
	byte *m_RevA;
	byte *m_RevB;
	float *m_Buffer1;
	float *m_Buffer2;

public:
	XDPMem()
		{
		Clear(true);
		}

	~XDPMem()
		{
		Clear(false);
		}

	void Clear(bool ctor = false)
		{
		if (!ctor)
			{
			m_TBBit.Clear();
			myfree(m_Buffer1);
			myfree(m_Buffer2);
			myfree(m_RevA);
			myfree(m_RevB);
			}
		
		m_MaxLA = 0;
		m_MaxLB = 0;
		m_RevA = 0;
		m_RevB = 0;
		m_Buffer1 = 0;
		m_Buffer2 = 0;
		}

	byte *GetRevA()
		{
		return m_RevA;
		}

	byte *GetRevB()
		{
		return m_RevB;
		}

	byte **GetTBBit()
		{
		return m_TBBit.GetData();
		}

	float *GetDPRow1()
		{
		return m_Buffer1 + 1;
		}

	float *GetDPRow2()
		{
		return m_Buffer2 + 1;
		}

	void Alloc(unsigned LA, unsigned LB)
		{
		if (LA + 8 > m_TBBit.m_AllocatedRowCount || LB + 8 > m_TBBit.m_AllocatedColCount)
			{
			m_TBBit.Alloc("TBBit", LA+129, LB+129);
			}

		if (LA > m_MaxLA)
			{
			m_MaxLA = LA + 128;
			myfree(m_RevA);
			m_RevA = myalloc(byte, m_MaxLA);
			}

		if (LB > m_MaxLB)
			{
			myfree(m_Buffer1);
			myfree(m_Buffer2);
			myfree(m_RevB);

			m_MaxLB = LB + 128;

		// Allow use of [-1]
			m_Buffer1 = myalloc(float, m_MaxLB+3);
			m_Buffer2 = myalloc(float, m_MaxLB+3);
			m_RevB = myalloc(byte, m_MaxLB);
			}
		asserta(m_TBBit.m_AllocatedRowCount >= LA && m_TBBit.m_AllocatedColCount >= LB);
		}

	void LogAlloc() const
		{
		Log("XDPMem[%p] m_MaxLA %u, m_MaxLB %u\n", this, m_MaxLA, m_MaxLB);
		Log("m_TBBit Rows %u, Cols %u, AllocRows %u, AllocCols %u\n",
		  m_TBBit.m_RowCount, m_TBBit.m_ColCount,
		  m_TBBit.m_AllocatedRowCount, m_TBBit.m_AllocatedColCount);
		}
	};

void LogAln(const byte *A, const byte *B, const char *Path);
class PathInfo;
void TraceBackBitMem(XDPMem &Mem, uint LA, uint LB, char State, PathInfo &PI);

#endif // xdpmem_h
