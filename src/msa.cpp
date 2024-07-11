#include "muscle.h"
#include "msa.h"
#include "textfile.h"
#include "seq.h"
#include "sequence.h"
#include <math.h>

const unsigned DEFAULT_SEQ_LENGTH = 500;

unsigned MSA::m_uIdCount = 0;

MSA::MSA()
	{
	m_uSeqCount = 0;
	m_uColCount = 0;

	m_szSeqs = 0;
	m_szNames = 0;

	m_IdToSeqIndex = 0;
	m_SeqIndexToId = 0;

	m_uCacheSeqCount = 0;
	m_uCacheSeqLength = 0;
	}

MSA::~MSA()
	{
	Free();
	}

void MSA::Free()
	{
	for (unsigned n = 0; n < m_uSeqCount; ++n)
		{
		delete[] m_szSeqs[n];
		delete[] m_szNames[n];
		}

	delete[] m_szSeqs;
	delete[] m_szNames;
	delete[] m_IdToSeqIndex;
	delete[] m_SeqIndexToId;

	m_uSeqCount = 0;
	m_uColCount = 0;

	m_szSeqs = 0;
	m_szNames = 0;

	m_IdToSeqIndex = 0;
	m_SeqIndexToId = 0;

	m_uCacheSeqLength = 0;
	m_uCacheSeqCount = 0;
	}

void MSA::SetSize(unsigned uSeqCount, unsigned uColCount)
	{
	Free();

	m_uSeqCount = uSeqCount;
	m_uCacheSeqLength = uColCount;
	m_uColCount = uColCount;

	if (0 == uSeqCount && 0 == uColCount)
		return;

	m_szSeqs = new char *[uSeqCount];
	m_szNames = new char *[uSeqCount];

	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		m_szSeqs[uSeqIndex] = new char[uColCount+1];
		m_szNames[uSeqIndex] = 0;
#if	DEBUG
		memset(m_szSeqs[uSeqIndex], '?', uColCount);
#endif
		m_szSeqs[uSeqIndex][uColCount] = 0;
		}

	if (m_uIdCount > 0)
		{
		m_IdToSeqIndex = new unsigned[m_uIdCount];
		m_SeqIndexToId = new unsigned[m_uSeqCount];
#if	DEBUG
		memset(m_IdToSeqIndex, 0xff, m_uIdCount*sizeof(unsigned));
		memset(m_SeqIndexToId, 0xff, m_uSeqCount*sizeof(unsigned));
#endif
		}
	}

void MSA::LogMe() const
	{
	if (0 == GetColCount())
		{
		Log("MSA empty\n");
		return;
		}

	const unsigned uColsPerLine = 50;
	unsigned uLinesPerSeq = (GetColCount() - 1)/uColsPerLine + 1;
	for (unsigned n = 0; n < uLinesPerSeq; ++n)
		{
		unsigned i;
		unsigned iStart = n*uColsPerLine;
		unsigned iEnd = GetColCount();
		if (iEnd - iStart + 1 > uColsPerLine)
			iEnd = iStart + uColsPerLine;
		Log("                       ");
		for (i = iStart; i < iEnd; ++i)
			Log("%u", i%10);
		Log("\n");
		Log("                       ");
		for (i = iStart; i + 9 < iEnd; i += 10)
			Log("%-10u", i);
		if (n == uLinesPerSeq - 1)
			Log(" %-10u", GetColCount());
		Log("\n");
		for (unsigned uSeqIndex = 0; uSeqIndex < m_uSeqCount; ++uSeqIndex)
			{
			Log("%12.12s", m_szNames[uSeqIndex]);
			Log("        ");
			Log("   ");
			for (i = iStart; i < iEnd; ++i)
				Log("%c", GetChar(uSeqIndex, i));
			if (0 != m_SeqIndexToId)
				Log(" [%5u]", m_SeqIndexToId[uSeqIndex]);
			Log("\n");
			}
		Log("\n\n");
		}
	}

char MSA::GetChar(unsigned uSeqIndex, unsigned uIndex) const
	{
// TODO: Performance cost?
	if (uSeqIndex >= m_uSeqCount || uIndex >= m_uColCount)
		Die("MSA::GetChar(%u/%u,%u/%u)",
		  uSeqIndex, m_uSeqCount, uIndex, m_uColCount);

	char c = m_szSeqs[uSeqIndex][uIndex];
//	assert(IsLegalChar(c));
	return c;
	}

unsigned MSA::GetLetter(unsigned uSeqIndex, unsigned uIndex) const
	{
// TODO: Performance cost?
	char c = GetChar(uSeqIndex, uIndex);
	unsigned uLetter = CharToLetter(c);
	if (uLetter >= 20)
		{
		char c = ' ';
		if (uSeqIndex < m_uSeqCount && uIndex < m_uColCount)
			c = m_szSeqs[uSeqIndex][uIndex];
		Die("MSA::GetLetter(%u/%u, %u/%u)='%c'/%u",
		  uSeqIndex, m_uSeqCount, uIndex, m_uColCount, c, uLetter);
		}
	return uLetter;
	}

//unsigned MSA::GetLetterEx(unsigned uSeqIndex, unsigned uIndex) const
//	{
//// TODO: Performance cost?
//	char c = GetChar(uSeqIndex, uIndex);
//	unsigned uLetter = CharToLetterEx(c);
//	return uLetter;
//	}

void MSA::SetSeqName(unsigned uSeqIndex, const char szName[])
	{
	if (uSeqIndex >= m_uSeqCount)
		Die("MSA::SetSeqName(%u, %s), count=%u", uSeqIndex, m_uSeqCount);
	delete[] m_szNames[uSeqIndex];
	int n = (int) strlen(szName) + 1;
	m_szNames[uSeqIndex] = new char[n];
	memcpy(m_szNames[uSeqIndex], szName, n);
	}

void MSA::GetSeqLabel(uint SeqIndex, string &Label) const
	{
	Label = string(GetSeqName(SeqIndex));
	}

const char *MSA::GetSeqName(unsigned uSeqIndex) const
	{
	if (uSeqIndex >= m_uSeqCount)
		Die("MSA::GetSeqName(%u), count=%u", uSeqIndex, m_uSeqCount);
	return m_szNames[uSeqIndex];
	}

bool MSA::IsGap(unsigned uSeqIndex, unsigned uIndex) const
	{
	char c = GetChar(uSeqIndex, uIndex);
	return IsGapChar(c);
	}

//bool MSA::IsWildcard(unsigned uSeqIndex, unsigned uIndex) const
//	{
//	char c = GetChar(uSeqIndex, uIndex);
//	return IsWildcardChar(c);
//	}

void MSA::SetChar(unsigned uSeqIndex, unsigned uIndex, char c)
	{
	if (uSeqIndex >= m_uSeqCount || uIndex > m_uCacheSeqLength)
		Die("MSA::SetChar(%u,%u)", uSeqIndex, uIndex);

	if (uIndex == m_uCacheSeqLength)
		{
		const unsigned uNewCacheSeqLength = m_uCacheSeqLength + DEFAULT_SEQ_LENGTH;
		for (unsigned n = 0; n < m_uSeqCount; ++n)
			{
			char *ptrNewSeq = new char[uNewCacheSeqLength+1];
			memcpy(ptrNewSeq, m_szSeqs[n], m_uCacheSeqLength);
			memset(ptrNewSeq + m_uCacheSeqLength, '?', DEFAULT_SEQ_LENGTH);
			ptrNewSeq[uNewCacheSeqLength] = 0;
			delete[] m_szSeqs[n];
			m_szSeqs[n] = ptrNewSeq;
			}

		m_uColCount = uIndex;
		m_uCacheSeqLength = uNewCacheSeqLength;
		}

	if (uIndex >= m_uColCount)
		m_uColCount = uIndex + 1;
	m_szSeqs[uSeqIndex][uIndex] = c;
	}

const char *MSA::GetSeqCharPtr(uint SeqIndex) const
	{
	asserta(SeqIndex < m_uSeqCount);
	const char *SeqCharPtr = m_szSeqs[SeqIndex];
	return SeqCharPtr;
	}

void MSA::GetRowStr(unsigned uSeqIndex, string &RowStr) const
	{
	RowStr.clear();
	const char *SeqCharPtr = GetSeqCharPtr(uSeqIndex);
	for (uint i = 0; i < m_uColCount; ++i)
		RowStr += SeqCharPtr[i];
	}

void MSA::GetSeq(unsigned uSeqIndex, Seq &seq) const
	{
	assert(uSeqIndex < m_uSeqCount);

	seq.Clear();

	for (unsigned n = 0; n < m_uColCount; ++n)
		if (!IsGap(uSeqIndex, n))
			{
			char c = GetChar(uSeqIndex, n);
			if (!isalpha(c))
				Die("Invalid character '%c' in sequence", c);
			c = toupper(c);
			seq.push_back(c);
			}
	const char *ptrName = GetSeqName(uSeqIndex);
	seq.SetName(ptrName);
	}

bool MSA::HasGap() const
	{
	for (unsigned uSeqIndex = 0; uSeqIndex < GetSeqCount(); ++uSeqIndex)
		for (unsigned n = 0; n < GetColCount(); ++n)
			if (IsGap(uSeqIndex, n))
				return true;
	return false;
	}

bool MSA::IsLegalLetter(unsigned uLetter) const
	{
	return uLetter < 20;
	}

void MSA::SetSeqCount(unsigned uSeqCount)
	{
	Free();
	SetSize(uSeqCount, DEFAULT_SEQ_LENGTH);
	}

void MSA::CopyCol(unsigned uFromCol, unsigned uToCol)
	{
	assert(uFromCol < GetColCount());
	assert(uToCol < GetColCount());
	if (uFromCol == uToCol)
		return;

	for (unsigned uSeqIndex = 0; uSeqIndex < GetSeqCount(); ++uSeqIndex)
		{
		const char c = GetChar(uSeqIndex, uFromCol);
		SetChar(uSeqIndex, uToCol, c);
		}
	}

void MSA::Copy(const MSA &msa)
	{
	Free();
	const unsigned uSeqCount = msa.GetSeqCount();
	const unsigned uColCount = msa.GetColCount();
	SetSize(uSeqCount, uColCount);

	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		SetSeqName(uSeqIndex, msa.GetSeqName(uSeqIndex));
		const unsigned uId = msa.GetSeqId(uSeqIndex);
		if (uId != UINT_MAX)
			SetSeqId(uSeqIndex, uId);
		for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
			{
			const char c = msa.GetChar(uSeqIndex, uColIndex);
			SetChar(uSeqIndex, uColIndex, c);
			}
		}
	}

void MSA::GetUpperLowerGapCount(uint ColIndex,
  uint &NU, uint &NL, uint &NG, uint &NDots, uint &NDashes) const
	{
	NU = 0;
	NL = 0;
	NG = 0;
	NDots = 0;
	NDashes = 0;
	for (unsigned SeqIndex = 0; SeqIndex < GetSeqCount(); ++SeqIndex)
		{
		char c = GetChar(SeqIndex, ColIndex);
		if (c == '.')
			{
			++NDots;
			++NG;
			}
		else if (c == '-')
			{
			++NDashes;
			++NG;
			}
		else if (isupper(c))
			++NU;
		else if (islower(c))
			++NL;
		}
	}

uint MSA::GetGapCount(unsigned uColIndex) const
	{
	uint n = 0;
	for (unsigned uSeqIndex = 0; uSeqIndex < GetSeqCount(); ++uSeqIndex)
		if (IsGap(uSeqIndex, uColIndex))
			++n;
	return n;
	}

bool MSA::IsGapColumn(unsigned uColIndex) const
	{
	assert(GetSeqCount() > 0);
	for (unsigned uSeqIndex = 0; uSeqIndex < GetSeqCount(); ++uSeqIndex)
		if (!IsGap(uSeqIndex, uColIndex))
			return false;
	return true;
	}

uint MSA::GetSeqIndex(const string &Label, bool FailOnError) const
	{
	for (unsigned uSeqIndex = 0; uSeqIndex < GetSeqCount(); ++uSeqIndex)
		if (0 == stricmp(Label.c_str(), GetSeqName(uSeqIndex)))
			return uSeqIndex;
	if (FailOnError)
		Die("Not found >%s", Label.c_str());
	return UINT_MAX;
	}

bool MSA::GetSeqIndex(const char *ptrSeqName, unsigned *ptruSeqIndex) const
	{
	for (unsigned uSeqIndex = 0; uSeqIndex < GetSeqCount(); ++uSeqIndex)
		if (0 == stricmp(ptrSeqName, GetSeqName(uSeqIndex)))
			{
			*ptruSeqIndex = uSeqIndex;
			return true;
			}
	return false;
	}

void MSA::DeleteCol(unsigned uColIndex)
	{
	assert(uColIndex < m_uColCount);
	size_t n = m_uColCount - uColIndex;
	if (n > 0)
		{
		for (unsigned uSeqIndex = 0; uSeqIndex < GetSeqCount(); ++uSeqIndex)
			{
			char *ptrSeq = m_szSeqs[uSeqIndex];
			memmove(ptrSeq + uColIndex, ptrSeq + uColIndex + 1, n);
			}
		}
	--m_uColCount;
	}

void MSA::DeleteColumns(unsigned uColIndex, unsigned uColCount)
	{
	for (unsigned n = 0; n < uColCount; ++n)
		DeleteCol(uColIndex);
	}

void MSA::FromFile(TextFile &File)
	{
	FromFASTAFile(File);
	}

static void FmtChar(char c, unsigned uWidth)
	{
	Log("%c", c);
	for (unsigned n = 0; n < uWidth - 1; ++n)
		Log(" ");
	}

static void FmtInt(unsigned u, unsigned uWidth)
	{
	static char szStr[1024];
	assert(uWidth < sizeof(szStr));
	if (u > 0)
		sprintf(szStr, "%u", u);
	else
		strcpy(szStr, ".");
	Log(szStr);
	unsigned n = (unsigned) strlen(szStr);
	if (n < uWidth)
		for (unsigned i = 0; i < uWidth - n; ++i)
			Log(" ");
	}

static void FmtInt0(unsigned u, unsigned uWidth)
	{
	static char szStr[1024];
	assert(uWidth < sizeof(szStr));
	sprintf(szStr, "%u", u);
	Log(szStr);
	unsigned n = (unsigned) strlen(szStr);
	if (n < uWidth)
		for (unsigned i = 0; i < uWidth - n; ++i)
			Log(" ");
	}

static void FmtPad(unsigned n)
	{
	for (unsigned i = 0; i < n; ++i)
		Log(" ");
	}

void MSA::GetLabelToSeqIndex(vector<string> &Labels,
  map<string, uint> &LabelToSeqIndex) const
	{
	Labels.clear();
	LabelToSeqIndex.clear();
	for (uint SeqIndex = 0; SeqIndex < m_uSeqCount; ++SeqIndex)
		{
		const string Label = (string) GetSeqName(SeqIndex);
		if (LabelToSeqIndex.find(Label) != LabelToSeqIndex.end())
			Die("Dupe  label >%s", Label.c_str());

		Labels.push_back(Label);
		LabelToSeqIndex[Label] = SeqIndex;
		}
	}

void MSA::FromSequence(const Sequence &s)
	{
	unsigned uSeqLength = s.GetLength();
	SetSize(1, uSeqLength);
	SetSeqName(0, s.GetLabelCStr());
	const byte *CharSeq = s.GetBytePtr();
	for (unsigned n = 0; n < uSeqLength; ++n)
		SetChar(0, n, CharSeq[n]);
	}

void MSA::FromMultiSequence(const MultiSequence &MS)
	{
	Clear();

	m_uSeqCount = MS.GetSeqCount();
	m_uColCount = MS.GetColCount();

	m_szSeqs = myalloc(char *, m_uSeqCount);
	m_szNames = myalloc(char *, m_uSeqCount);
	for (uint i = 0; i < m_uSeqCount; ++i)
		{
		m_szSeqs[i] = (char *) MS.GetCharPtr(i);
		m_szNames[i] = (char *) MS.GetLabel(i);
		}
	}

void MSA::FromSeq(const Seq &s)
	{
	unsigned uSeqLength = s.Length();
	SetSize(1, uSeqLength);
	SetSeqName(0, s.GetName());
	if (0 != m_SeqIndexToId)
		SetSeqId(0, s.GetId());
	for (unsigned n = 0; n < uSeqLength; ++n)
		SetChar(0, n, s[n]);
	}

unsigned MSA::GetCharCount(unsigned uSeqIndex, unsigned uColIndex) const
	{
	assert(uSeqIndex < GetSeqCount());
	assert(uColIndex < GetColCount());

	unsigned uCol = 0;
	for (unsigned n = 0; n <= uColIndex; ++n)
		if (!IsGap(uSeqIndex, n))
			++uCol;
	return uCol;
	}

void MSA::CopySeq(unsigned uToSeqIndex, const MSA &msaFrom, unsigned uFromSeqIndex)
	{
	assert(uToSeqIndex < m_uSeqCount);
	const unsigned uColCount = msaFrom.GetColCount();
	assert(m_uColCount == uColCount ||
	  (0 == m_uColCount && uColCount <= m_uCacheSeqLength));

	memcpy(m_szSeqs[uToSeqIndex], msaFrom.GetSeqBuffer(uFromSeqIndex), uColCount);
	SetSeqName(uToSeqIndex, msaFrom.GetSeqName(uFromSeqIndex));
	if (0 == m_uColCount)
		m_uColCount = uColCount;
	}

const char *MSA::GetSeqBuffer(unsigned uSeqIndex) const
	{
	assert(uSeqIndex < m_uSeqCount);
	return m_szSeqs[uSeqIndex];
	}

void MSA::DeleteSeq(unsigned uSeqIndex)
	{
	assert(uSeqIndex < m_uSeqCount);

	delete m_szSeqs[uSeqIndex];
	delete m_szNames[uSeqIndex];

	const unsigned uBytesToMove = (m_uSeqCount - uSeqIndex)*sizeof(char *);
	if (uBytesToMove > 0)
		{
		memmove(m_szSeqs + uSeqIndex, m_szSeqs + uSeqIndex + 1, uBytesToMove);
		memmove(m_szNames + uSeqIndex, m_szNames + uSeqIndex + 1, uBytesToMove);
		}

	--m_uSeqCount;
	}

bool MSA::IsEmptyCol(unsigned uColIndex) const
	{
	const unsigned uSeqCount = GetSeqCount();
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		if (!IsGap(uSeqIndex, uColIndex))
			return false;
	return true;
	}

//void MSA::DeleteEmptyCols(bool bProgress)
//	{
//	unsigned uColCount = GetColCount();
//	for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
//		{
//		if (IsEmptyCol(uColIndex))
//			{
//			if (bProgress)
//				{
//				Log("Deleting col %u of %u\n", uColIndex, uColCount);
//				printf("Deleting col %u of %u\n", uColIndex, uColCount);
//				}
//			DeleteCol(uColIndex);
//			--uColCount;
//			}
//		}
//	}

unsigned MSA::AlignedColIndexToColIndex(unsigned uAlignedColIndex) const
	{
	Die("MSA::AlignedColIndexToColIndex not implemented");
	return 0;
	}

bool MSA::SeqsEq(const MSA &a1, unsigned uSeqIndex1, const MSA &a2,
  unsigned uSeqIndex2)
	{
	Seq s1;
	Seq s2;

	a1.GetSeq(uSeqIndex1, s1);
	a2.GetSeq(uSeqIndex2, s2);

	s1.StripGaps();
	s2.StripGaps();

	return s1.EqIgnoreCase(s2);
	}

unsigned MSA::GetUngappedSeqLength(unsigned uSeqIndex) const
	{
	assert(uSeqIndex < GetSeqCount());

	const unsigned uColCount = GetColCount();
	unsigned uLength = 0;
	for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
		if (!IsGap(uSeqIndex, uColIndex))
			++uLength;
	return uLength;
	}

void MSA::GetPWID(unsigned uSeqIndex1, unsigned uSeqIndex2, double *ptrPWID,
  unsigned *ptruPosCount) const
	{
	assert(uSeqIndex1 < GetSeqCount());
	assert(uSeqIndex2 < GetSeqCount());

	unsigned uSameCount = 0;
	unsigned uPosCount = 0;
	const unsigned uColCount = GetColCount();
	for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
		{
		char c1 = GetChar(uSeqIndex1, uColIndex);
		if (IsGapChar(c1))
			continue;
		char c2 = GetChar(uSeqIndex2, uColIndex);
		if (IsGapChar(c2))
			continue;
		++uPosCount;
		if (c1 == c2)
			++uSameCount;
		}
	*ptruPosCount = uPosCount;
	if (uPosCount > 0)
		*ptrPWID = 100.0 * (double) uSameCount / (double) uPosCount;
	else
		*ptrPWID = 0;
	}

//unsigned MSA::UniqueResidueTypes(unsigned uColIndex) const
//	{
//	assert(uColIndex < GetColCount());
//
//	unsigned Counts[MAX_ALPHA];
//	memset(Counts, 0, sizeof(Counts));
//	const unsigned uSeqCount = GetSeqCount();
//	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
//		{
//		if (IsGap(uSeqIndex, uColIndex) || IsWildcard(uSeqIndex, uColIndex))
//			continue;
//		const unsigned uLetter = GetLetter(uSeqIndex, uColIndex);
//		++(Counts[uLetter]);
//		}
//	unsigned uUniqueCount = 0;
//	for (unsigned uLetter = 0; uLetter < g_AlphaSize; ++uLetter)
//		if (Counts[uLetter] > 0)
//			++uUniqueCount;
//	return uUniqueCount;
//	}

double MSA::GetOcc(unsigned uColIndex) const
	{
	unsigned uGapCount = 0;
	for (unsigned uSeqIndex = 0; uSeqIndex < GetSeqCount(); ++uSeqIndex)
		if (IsGap(uSeqIndex, uColIndex))
			++uGapCount;
	unsigned uSeqCount = GetSeqCount();
	return (double) (uSeqCount - uGapCount) / (double) uSeqCount;
	}

void MSA::ToFile(TextFile &File) const
	{
	ToFASTAFile(File);
	}

void MSA::GetUngappedSeqStr(uint SeqIndex, string &SeqStr) const
	{
	SeqStr.clear();
	asserta(SeqIndex < m_uSeqCount);
	for (uint i = 0; i < m_uColCount; ++i)
		{
		char c = m_szSeqs[SeqIndex][i];
		if (!isgap(c))
			SeqStr += toupper(c);
		}
	}

bool MSA::ColumnHasGap(unsigned uColIndex) const
	{
	const unsigned uSeqCount = GetSeqCount();
	for (unsigned uSeqIndex = 0; uSeqIndex < uSeqCount; ++uSeqIndex)
		if (IsGap(uSeqIndex, uColIndex))
			return true;
	return false;
	}

void MSA::SetIdCount(unsigned uIdCount)
	{
	//if (m_uIdCount != 0)
	//	Die("MSA::SetIdCount: may only be called once");

	if (m_uIdCount > 0)
		{
		if (uIdCount > m_uIdCount)
			Die("MSA::SetIdCount: cannot increase count");
		return;
		}
	m_uIdCount = uIdCount;
	}

void MSA::SetSeqId(unsigned uSeqIndex, unsigned uId)
	{
	assert(uSeqIndex < m_uSeqCount);
	assert(uId == UINT_MAX || uId < m_uIdCount);
	if (0 == m_SeqIndexToId)
		{
		if (0 == m_uIdCount)
			Die("MSA::SetSeqId, SetIdCount has not been called");
		m_IdToSeqIndex = new unsigned[m_uIdCount];
		m_SeqIndexToId = new unsigned[m_uSeqCount];

		memset(m_IdToSeqIndex, 0xff, m_uIdCount*sizeof(unsigned));
		memset(m_SeqIndexToId, 0xff, m_uSeqCount*sizeof(unsigned));
		}
	m_SeqIndexToId[uSeqIndex] = uId;
	if (uId != UINT_MAX)
		m_IdToSeqIndex[uId] = uSeqIndex;
	}

unsigned MSA::GetSeqIndex(unsigned uId) const
	{
	assert(uId < m_uIdCount);
	assert(0 != m_IdToSeqIndex);
	unsigned uSeqIndex = m_IdToSeqIndex[uId];
	assert(uSeqIndex < m_uSeqCount);
	return uSeqIndex;
	}

bool MSA::GetSeqIndex(unsigned uId, unsigned *ptruIndex) const
	{
	for (unsigned uSeqIndex = 0; uSeqIndex < m_uSeqCount; ++uSeqIndex)
		{
		if (uId == m_SeqIndexToId[uSeqIndex])
			{
			*ptruIndex = uSeqIndex;
			return true;
			}
		}
	return false;
	}

unsigned MSA::GetSeqId(unsigned uSeqIndex) const
	{
	if (m_SeqIndexToId == 0)
		return UINT_MAX;
	assert(uSeqIndex < m_uSeqCount);
	unsigned uId = m_SeqIndexToId[uSeqIndex];
	assert(uId == UINT_MAX || uId < m_uIdCount);
	return uId;
	}

void MSASubsetByIds(const MSA &msaIn, const unsigned Ids[], unsigned uIdCount,
  MSA &msaOut)
	{
	const unsigned uColCount = msaIn.GetColCount();
	msaOut.SetSize(uIdCount, uColCount);
	for (unsigned uSeqIndexOut = 0; uSeqIndexOut < uIdCount; ++uSeqIndexOut)
		{
		const unsigned uId = Ids[uSeqIndexOut];

		const unsigned uSeqIndexIn = msaIn.GetSeqIndex(uId);
		const char *ptrName = msaIn.GetSeqName(uSeqIndexIn);

		msaOut.SetSeqId(uSeqIndexOut, uId);
		msaOut.SetSeqName(uSeqIndexOut, ptrName);

		for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
			{
			const char c = msaIn.GetChar(uSeqIndexIn, uColIndex);
			msaOut.SetChar(uSeqIndexOut, uColIndex, c);
			}
		}
	}

// Caller must allocate ptrSeq and ptrLabel as new char[n].
void MSA::AppendSeq(char *ptrSeq, unsigned uSeqLength, char *ptrLabel)
	{
	if (m_uSeqCount > m_uCacheSeqCount)
		Die("Internal error MSA::AppendSeq");
	if (m_uSeqCount == m_uCacheSeqCount)
		ExpandCache(m_uSeqCount + 4, uSeqLength);
	m_szSeqs[m_uSeqCount] = ptrSeq;
	m_szNames[m_uSeqCount] = ptrLabel;
	++m_uSeqCount;
	}

void MSA::ExpandCache(unsigned uSeqCount, unsigned uColCount)
	{
	if (m_IdToSeqIndex != 0 || m_SeqIndexToId != 0 || uSeqCount < m_uSeqCount)
		Die("Internal error MSA::ExpandCache");

	if (m_uSeqCount > 0 && uColCount != m_uColCount)
		Die("Internal error MSA::ExpandCache, ColCount changed");

	char **NewSeqs = new char *[uSeqCount];
	char **NewNames = new char *[uSeqCount];

	for (unsigned uSeqIndex = 0; uSeqIndex < m_uSeqCount; ++uSeqIndex)
		{
		NewSeqs[uSeqIndex] = m_szSeqs[uSeqIndex];
		NewNames[uSeqIndex] = m_szNames[uSeqIndex];
		}

	for (unsigned uSeqIndex = m_uSeqCount; uSeqIndex < uSeqCount; ++uSeqIndex)
		{
		char *Seq = new char[uColCount];
		NewSeqs[uSeqIndex] = Seq;
#if	DEBUG
		memset(Seq, '?', uColCount);
#endif
		}

	delete[] m_szSeqs;
	delete[] m_szNames;

	m_szSeqs = NewSeqs;
	m_szNames = NewNames;

	m_uCacheSeqCount = uSeqCount;
	m_uCacheSeqLength = uColCount;
	m_uColCount = uColCount;
	}

//void MSA::FixAlpha()
//	{
//	ClearInvalidLetterWarning();
//	for (unsigned uSeqIndex = 0; uSeqIndex < m_uSeqCount; ++uSeqIndex)
//		{
//		for (unsigned uColIndex = 0; uColIndex < m_uColCount; ++uColIndex)
//			{
//			char c = GetChar(uSeqIndex, uColIndex);
//			if (!IsResidueChar(c) && !IsGapChar(c))
//				{
//				char w = GetWildcardChar();
//				// Warning("Invalid letter '%c', replaced by '%c'", c, w);
//				InvalidLetterWarning(c, w);
//				SetChar(uSeqIndex, uColIndex, w);
//				}
//			}
//		}
//	ReportInvalidLetters();
//	}

//ALPHA MSA::GuessAlpha() const
//	{
//// If at least MIN_NUCLEO_PCT of the first CHAR_COUNT non-gap
//// letters belong to the nucleotide alphabet, guess nucleo.
//// Otherwise amino.
//	const unsigned CHAR_COUNT = 100;
//	const unsigned MIN_NUCLEO_PCT = 95;
//
//	const unsigned uSeqCount = GetSeqCount();
//	const unsigned uColCount = GetColCount();
//	if (0 == uSeqCount)
//		return ALPHA_Amino;
//
//	unsigned uNucleoCount = 0;
//	unsigned uTotal = 0;
//	unsigned i = 0;
//	for (;;)
//		{
//		unsigned uSeqIndex = i/uColCount;
//		if (uSeqIndex >= uSeqCount)
//			break;
//		unsigned uColIndex = i%uColCount;
//		++i;
//		char c = GetChar(uSeqIndex, uColIndex);
//		if (IsGapChar(c))
//			continue;
//		if (g_CharToLetterNucleo[c] < 4 || c == 'N' || c == 'n')
//			++uNucleoCount;
//		++uTotal;
//		if (uTotal >= CHAR_COUNT)
//			break;
//		}
//	if (uTotal != 0 && ((uNucleoCount*100)/uTotal) >= MIN_NUCLEO_PCT)
//		return ALPHA_Nucleo;
//	return ALPHA_Amino;
//	}

void MSA::GetPosToCol(uint SeqIndex, vector<uint> &PosToCol) const
	{
	PosToCol.clear();
	const uint ColCount = GetColCount();
	const char *Seq = GetSeqCharPtr(SeqIndex);
	PosToCol.reserve(ColCount);
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = Seq[Col];
		if (!isgap(c))
			PosToCol.push_back(Col);
		}
	}

// 1-based positions, if <0 the column has a gap in this
// sequence which follows at 1-based position (-Pos).
// If ColToPos[Col] is 0, this is left terminal gap.
void MSA::GetColToPos1(uint SeqIndex, vector<int> &ColToPos) const
	{
	ColToPos.clear();
	const uint ColCount = GetColCount();
	const char *Seq = GetSeqCharPtr(SeqIndex);
	ColToPos.reserve(ColCount);
	int Pos = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = Seq[Col];
		if (isgap(c))
			ColToPos.push_back(Pos == 0 ? 0 : -Pos);
		else
			ColToPos.push_back(++Pos);
		}
	}

void MSA::GetColToPos(uint SeqIndex, vector<uint> &ColToPos) const
	{
	ColToPos.clear();
	const uint ColCount = GetColCount();
	const char *Seq = GetSeqCharPtr(SeqIndex);
	ColToPos.reserve(ColCount);
	uint Pos = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = Seq[Col];
		if (isgap(c))
			ColToPos.push_back(UINT_MAX);
		else
			ColToPos.push_back(Pos++);
		}
	}

bool MSA::ColIsUpper(uint ColIndex, double MaxGapFract) const
	{
	const uint SeqCount = GetSeqCount();
	uint UpperCount = 0;
	uint LowerCount = 0;
	uint GapCount = 0;
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		char c = GetChar(SeqIndex, ColIndex);
		if (isgap(c))
			{
			++GapCount;
			continue;
			}
		if (!isalpha(c))
			continue;
		if (isupper(c))
			++UpperCount;
		else
			++LowerCount;
		}

	if (UpperCount == 0 && LowerCount == 0)
		return false;

	if (UpperCount > 0 && LowerCount > 0)
		Die("Column %u has mixed case letters", ColIndex);

	if (double(GapCount)/SeqCount > MaxGapFract)
		return false;

	if (UpperCount == 0)
		return false;

	return true;
	}

bool MSA::ColIsAligned(uint ColIndex) const
	{
	const uint SeqCount = GetSeqCount();
	uint UpperCount = 0;
	uint LowerCount = 0;
	uint GapCount = 0;
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		char c = GetChar(SeqIndex, ColIndex);
		if (isgap(c))
			{
			++GapCount;
			continue;
			}
		if (!isalpha(c))
			continue;
		if (isupper(c))
			++UpperCount;
		else
			++LowerCount;
		}

	if (UpperCount == 0 && LowerCount == 0)
		return false;

	if (UpperCount > 0 && LowerCount > 0)
		Die("Column %u has mixed case letters", ColIndex);

	if (UpperCount < 2)
		return false;

	return true;
	}

void MSA::DeleteAllGapCols(MSA &msa) const
	{
	const uint ColCount = GetColCount();
	const uint SeqCount = GetSeqCount();
	uint NewColCount = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		uint GapCount = GetGapCount(Col);
		if (GapCount != SeqCount)
			++NewColCount;
		}
	msa.SetSize(SeqCount, NewColCount);
	uint NewColIndex = 0;
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		msa.m_szNames[SeqIndex] = mystrsave(GetLabel(SeqIndex));

	for (uint Col = 0; Col < ColCount; ++Col)
		{
		uint GapCount = GetGapCount(Col);
		if (GapCount == SeqCount)
			continue;
		for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
			{
			char c = GetChar(SeqIndex, Col);
			msa.SetChar(SeqIndex, NewColIndex, c);
			}
		++NewColIndex;
		}
	}
