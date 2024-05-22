#ifndef pathinfo_h
#define pathinfo_h

#include "obj.h"
#include "xdpmem.h"

class PathInfo : public Obj
	{
	friend class ObjMgr;

private:
	unsigned m_BufferBytes;
	unsigned m_ColCount;
	char *m_Buffer;

protected:
	PathInfo() : Obj(OT_PathInfo)
		{
		m_BufferBytes = 0;
		m_Buffer = 0;
		}

	virtual ~PathInfo()
		{
		myfree(m_Buffer);
		m_Buffer = 0;
		}

public:
	virtual void Reset()
		{
		SetEmpty();
		}

	virtual unsigned GetMemBytes() const
		{
		return m_BufferBytes;
		}

public:
	void Realloc(unsigned Bytes);
	void SetEmpty();
	void FreeIfBig();
	void Alloc(unsigned Bytes);
	void Alloc2(unsigned LA, unsigned LB);
	void AppendPath(const PathInfo &PI);
	void PrependPath(const PathInfo &PI);
	void AppendMs(unsigned Count);
	void AppendDs(unsigned Count);
	void AppendIs(unsigned Count);
	void ToOps(string &Ops, vector<uint> &Lengths) const;

	void Grow(unsigned n)
		{
		if (m_ColCount + n < m_BufferBytes)
			return;
		Realloc(RoundUp(m_ColCount + n, 4096));
		}

	unsigned GetColCount() const
		{
		assert(m_ColCount == strlen(m_Buffer));
		return m_ColCount;
		}

	const char *GetPath() const
		{
		asserta(m_Buffer != 0);
		return m_Buffer;
		}

	void GetPathStr(string &PathStr) const
		{
		PathStr.clear();
		for (uint i = 0; i < m_ColCount; ++i)
			{
			char c = m_Buffer[i];
			asserta(c == 'M' || c == 'D' || c == 'I');
			PathStr += c;
			}
		}

	unsigned GetLeftICount();
	unsigned TrimLeftIs();
	void TrimRightIs();

	unsigned GetCounts(unsigned &M, unsigned &D, unsigned &I) const;
	void Reverse();
	void AppendChar(char c);
	};

#endif // pathinfo_h
