#ifndef mx_h
#define mx_h

#include <list>
#include <limits.h>
#include <math.h>
#include "myutils.h"

const int OPT_LOG = 0x01;
const int OPT_EXP = 0x02;
const int OPT_ZERO_BASED = 0x04;
const float MINUS_INFINITY = -9e9f;
const float UNINIT = -8e8f;

struct SeqData;

template<class T> const char *TypeToStr(T t)
	{
	Die("Unspecialised TypeToStr() called");
	return 0;
	}

template<> inline const char *TypeToStr<unsigned short>(unsigned short f)
	{
	static char s[16];

	sprintf(s, "%12u", f);
	return s;
	}

template<> inline const char *TypeToStr<unsigned>(unsigned f)
	{
	static char s[16];

	sprintf(s, "%12u", f);
	return s;
	}

template<> inline const char *TypeToStr<short>(short f)
	{
	static char s[16];

	sprintf(s, "%12d", f);
	return s;
	}

template<> inline const char *TypeToStr<int>(int f)
	{
	static char s[16];

	sprintf(s, "%5d", f);
	return s;
	}

template<> inline const char *TypeToStr<float>(float f)
	{
	static char s[16];

	if (f == UNINIT)
		sprintf(s, "%12.12s", "?");
	else if (f < MINUS_INFINITY/2)
		sprintf(s, "%12.12s", "*");
	else if (f == 0.0f)
		sprintf(s, "%12.12s", ".");
	else if (f >= -1e5 && f <= 1e5)
		sprintf(s, "%12.5f", f);
	else
		sprintf(s, "%12.4g", f);
	return s;
	}

template<> inline const char *TypeToStr<double>(double f)
	{
	static char s[16];

	if (f < -1e9)
		sprintf(s, "%12.12s", "*");
	else if (f == 0.0f)
		sprintf(s, "%12.12s", ".");
	else if (f >= -1e-5 && f <= 1e5)
		sprintf(s, "%12.5f", f);
	else
		sprintf(s, "%12.4g", f);
	return s;
	}

static inline const char *FloatToStr(float f, string &s)
	{
	s = TypeToStr<float>(f);
	return s.c_str();
	}

template<> inline const char *TypeToStr<char>(char c)
	{
	static char s[2];
	s[0] = c;
	return s;
	}

template<> inline const char *TypeToStr<byte>(byte c)
	{
	static char s[2];
	s[0] = c;
	return s;
	}

template<> inline const char *TypeToStr<bool>(bool tof)
	{
	static char s[2];
	s[0] = tof ? 'T' : 'F';
	return s;
	}

struct MxBase
	{
private:
	MxBase(const MxBase &rhs);
	MxBase &operator=(const MxBase &rhs);

public:
	string m_Name;
	unsigned m_RowCount;
	unsigned m_ColCount;
	unsigned m_AllocatedRowCount;
	unsigned m_AllocatedColCount;

	static list<MxBase *> *m_Matrices;

	static unsigned m_AllocCount;
	static unsigned m_ZeroAllocCount;
	static unsigned m_GrowAllocCount;
	static double m_TotalBytes;
	static double m_MaxBytes;

	static void OnCtor(MxBase *Mx);
	static void OnDtor(MxBase *Mx);

	MxBase()
		{
		m_AllocatedRowCount = 0;
		m_AllocatedColCount = 0;
		m_RowCount = 0;
		m_ColCount = 0;
		OnCtor(this);
		}

	virtual ~MxBase()
		{
		OnDtor(this);
		}

public:
	virtual unsigned GetTypeSize() const = 0;
	virtual unsigned GetTotalBytes() const = 0;
	virtual void AllocData(unsigned RowCount, unsigned ColCount) = 0;
	virtual void FreeData() = 0;
	virtual const char *GetAsStr(unsigned i, unsigned j) const = 0;

public:
	void Clear()
		{
		FreeData();
		m_AllocatedRowCount = 0;
		m_AllocatedColCount = 0;
		m_RowCount = 0;
		m_ColCount = 0;
		}

	bool IsEmpty() const
		{
		return m_RowCount == 0;
		}

public:
	void Alloc(unsigned RowCount, unsigned ColCount, const string &Name = "");
	void Alloc(const string &Name, unsigned RowCount, unsigned ColCount)
		{
		Alloc(RowCount, ColCount, Name);
		}

	unsigned GetRowCount() const { return m_RowCount; }
	unsigned GetColCount() const { return m_ColCount; }

	static void LogAll()
		{
		Log("\n");
		if (m_Matrices == 0)
			{
			Log("MxBase::m_Matrices=0\n");
			return;
			}
		Log("\n");
		Log("AllRows  AllCols    Sz        MB  Name\n");
		Log("-------  -------  ----  --------  ----\n");
		double TotalMB = 0;
		for (list<MxBase *>::const_iterator p = m_Matrices->begin();
		  p != m_Matrices->end(); ++p)
			{
			const MxBase *Mx = *p;
			if (Mx == 0)
				continue;
			//if (Mx->m_RowCount != 0 || ShowEmpty)
			//	Mx->LogMe(WithData);
			unsigned ar = Mx->m_AllocatedRowCount;
			if (ar == 0)
				continue;
			unsigned ac = Mx->m_AllocatedColCount;
			unsigned sz = Mx->GetTypeSize();
			double MB = (double) ar*(double) ac*(double) sz/1e6;
			TotalMB += MB;
			Log("%7u  %7u  %4u  %8.2f  %s\n", ar, ac, sz, MB, Mx->m_Name.c_str());
			}
		Log("                        --------\n");
		Log("%7.7s  %7.7s  %4.4s  %8.2f\n", "", "", "", TotalMB);
		}

	void LogMe(bool WithData = true, int Opts = 0) const;
	static void LogCounts();
	};

template<class T> struct Mx : public MxBase
	{
// Disable unimplemented stuff
private:
	Mx(Mx &rhs);
	Mx &operator=(Mx &rhs);
	// const Mx &operator=(const Mx &rhs) const;

public:
	T **m_Data;

	Mx()
		{
		m_Data = 0;
		}
	
	~Mx()
		{
		FreeData();
		}

	virtual void AllocData(unsigned RowCount, unsigned ColCount)
		{
		double dRowPtrsBytes = double(sizeof(T *))*double(RowCount);
		if (dRowPtrsBytes > (double) (UINT_MAX-16))
			Die("Mx::AllocData Rows=%u, sizeof(T *)=%u, row %s bytes too big",
			  RowCount, sizeof(T *), MemBytesToStr(dRowPtrsBytes));

		double dDataBytes = double(sizeof(T))*double(RowCount)*double(ColCount);
		if (dDataBytes > (double) (UINT_MAX-16))
			Die("Mx::AllocData Rows=%u, Cols=%u, sizeof(T)=%u, data %s too big",
			  RowCount, ColCount, sizeof(T), MemBytesToStr(dDataBytes));

		unsigned RowBytes = sizeof(T *)*RowCount;
		unsigned DataBytes = sizeof(T)*RowCount*ColCount;

		byte *Buffer = myalloc(byte, RowBytes + DataBytes);
		m_Data = (T **) Buffer;
		byte *Base = Buffer + RowCount*sizeof(T *);
		for (unsigned i = 0; i < RowCount; ++i)
			m_Data[i] = (T *) (Base + i*ColCount*sizeof(T));

		//m_Data = myalloc(T *, RowCount);
		//for (unsigned i = 0; i < RowCount; ++i)
		//	m_Data[i] = myalloc(T, ColCount);

		m_AllocatedRowCount = RowCount;
		m_AllocatedColCount = ColCount;
		}

	virtual void FreeData()
		{
		//for (unsigned i = 0; i < m_AllocatedRowCount; ++i)
		//	myfree(m_Data[i]);
		//myfree(m_Data);
		myfree(m_Data);

		m_Data = 0;
		m_RowCount = 0;
		m_ColCount = 0;
		m_AllocatedRowCount = 0;
		m_AllocatedColCount = 0;
		}

	T **GetData()
		{
		return (T **) m_Data;
		}

	T Get(unsigned i, unsigned j) const
		{
		assert(i < m_RowCount);
		assert(j < m_ColCount);
		return m_Data[i][j];
		}

	void Put(unsigned i, unsigned j, T x)
		{
		assert(i < m_RowCount);
		assert(j < m_ColCount);
		m_Data[i][j] = x;
		}

	unsigned GetTypeSize() const
		{
		return sizeof(T);
		}

	virtual unsigned GetTotalBytes() const
		{
		return m_AllocatedRowCount*m_AllocatedColCount*GetTypeSize() +
		  m_AllocatedRowCount*sizeof(T *);
		}

	const char *GetAsStr(unsigned i, unsigned j) const
		{
		return TypeToStr<T>(Get(i, j));
		}

	const T *const *const GetData() const
		{
		return (const T *const *) m_Data;
		}

	void Zero()
		{
		PutAll(0);
		}

	void PutAll(T v)
		{
		for (unsigned i = 0; i < m_RowCount; ++i)
			for (unsigned j = 0; j < m_ColCount; ++j)
				m_Data[i][j] = v;
		}
	};

#endif // mx_h
