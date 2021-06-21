#ifndef diaglist_h
#define diaglist_h

const unsigned EMPTY = (unsigned) ~0;
const unsigned MAX_DIAGS = 1024;

struct Diag
	{
	unsigned m_uStartPosA;
	unsigned m_uStartPosB;
	unsigned m_uLength;
	};

struct Rect
	{
	unsigned m_uStartPosA;
	unsigned m_uStartPosB;
	unsigned m_uLengthA;
	unsigned m_uLengthB;
	};

class DiagList
	{
public:
	DiagList()
		{
		m_uCount = 0;
		}
	~DiagList()
		{
		Free();
		}

public:
// Creation
	void Clear()
		{
		Free();
		}
	void FromPath(const PWPath &Path);
	void Add(const Diag &d);
	void Add(unsigned uStartPosA, unsigned uStartPosB, unsigned uLength);
	void DeleteIncompatible();

// Accessors
	unsigned GetCount() const
		{
		return m_uCount;
		}
	const Diag &Get(unsigned uIndex) const;

// Operations
	void Sort();
	void Copy(const DiagList &DL);

// Query
	// returns true iff given diagonal is included in the list
	// in whole or in part.
	bool NonZeroIntersection(const Diag &d) const;
	bool IsSorted() const;

// Diagnostics
	void LogMe() const;

private:
	void Free()
		{
		m_uCount = 0;
		}

private:
	unsigned m_uCount;
	Diag m_Diags[MAX_DIAGS];
	};

unsigned DiagOverlap(const Diag &d1, const Diag &d2);
unsigned DiagOverlapA(const Diag &d1, const Diag &d2);
unsigned DiagOverlapB(const Diag &d1, const Diag &d2);
unsigned DiagBreak(const Diag &d1, const Diag &d2);
bool DiagCompatible(const Diag &d1, const Diag &d2);
void CheckDiags(const ProfPos *PA, unsigned uLengthA, const ProfPos *PB,
  unsigned uLengthB, const MSA &msaA, const MSA &msaB, const PWPath &Path);
void FindDiags(const ProfPos *PX, unsigned uLengthX, const ProfPos *PY,
  unsigned uLengthY, DiagList &DL);
void FindDiagsNuc(const ProfPos *PX, unsigned uLengthX, const ProfPos *PY,
  unsigned uLengthY, DiagList &DL);
void MergeDiags(DiagList &DL);

#endif // diaglist_h
