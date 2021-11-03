#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <cctype>
#include <cstdlib>
#include "filebuffer.h"

#define SEQ_TRACE	0

class Sequence;

// m_CharVec uses one-based indexing, m_CharVec[0] set to '@'.
class Sequence
	{
public:
#if SEQ_TRACE
	uint m_Id = UINT_MAX;
#endif

// Global input MulitSquence index
	uint m_GSI = UINT_MAX;

// Sparse matrix index
	uint m_SMI = UINT_MAX;

	string m_Label;
	vector<char> m_CharVec;

private:
	Sequence()
		{
#if SEQ_TRACE
		m_Id = UINT_MAX;
#endif
		m_Label = "~";
		m_GSI = UINT_MAX;
		m_SMI = UINT_MAX;
		}
	~Sequence() {}

public:
	bool FromFileBuffer(FileBuffer& infile, bool stripGaps = false);

	void Create(const vector<char> *m_CharVec, string m_Label,
	  uint GSI, uint SMI);

	void FromString(const string &Label, const string &Seq);
//	void Copy(const Sequence &rhs);

	void InitData()
		{
		m_CharVec.clear();
		m_CharVec.push_back('@');
		}

	void AppendChar(char c)
		{
		m_CharVec.push_back(c);
		}

	const string &GetLabel() const
		{
		return m_Label;
		}

	const char *GetLabelCStr() const
		{
		return m_Label.c_str();
		}

// One-based
	char *GetCharPtr1()
		{
		return m_CharVec.data();
		}

	const char *GetCharPtr1() const
		{
		return m_CharVec.data();
		}

	const byte *GetBytePtr() const
		{
		const char *CharPtr = m_CharVec.data() + 1;
		const byte *BytePtr = (const byte *) CharPtr;
		return BytePtr;
		}

// Chars stored with one-based indexing.
	char GetPosition(int i) const
		{
		assert(i >= 1 && i < m_CharVec.size());
		return m_CharVec[i];
		}

	char GetChar(uint ZeroBasedPos) const
		{
		return GetPosition(int(ZeroBasedPos+1));
		}

	void SetGSI(uint GSI)
		{
		asserta(m_GSI == UINT_MAX);
		m_GSI = GSI;
		}

	void OverwriteGSI(uint GSI)
		{
		m_GSI = GSI;
		}

	void OverwriteLabel(const string &Label)
		{
		m_Label = Label;
		}

	uint GetSMI() const
		{
		return m_SMI;
		}

	uint GetGSI() const
		{
		return m_GSI;
		}

	uint GetLength() const
		{
		uint n = (uint) m_CharVec.size();
		asserta(n > 0);
		uint L = n - 1;
		return L;
		}

	void WriteMFA(FILE *f) const;
	Sequence* Clone() const;
	Sequence* AddGaps(const vector<char>* alignment, char id) const;
	Sequence* AddGapsPath(const string &Path, char id) const;
	Sequence* DeleteGaps() const;
	void GetPosToCol_OneBased(vector<uint> &PosToCol) const;
	void GetPosToCol(vector<uint> &PosToCol) const;
	void GetColToPos(vector<uint> &ColToPos) const;
	void LogMe() const;
#if SEQ_TRACE
	void AssertId() const;
#endif

public:
#if SEQ_TRACE

	static Sequence *_NewSequence(const char *File, int Line);
#define NewSequence()	Sequence::_NewSequence(__FILE__, __LINE__)
	static void _DeleteSequence(const Sequence *Seq,
	  const char *File, int Line);
#define DeleteSequence(s)	Sequence::_DeleteSequence((s), __FILE__, __LINE__)
	static void AllocReport(const string &Msg);

#else

	static Sequence *_NewSequence();
#define NewSequence()	Sequence::_NewSequence()
	static void _DeleteSequence(const Sequence *Seq);
#define DeleteSequence(s)	Sequence::_DeleteSequence((s))

#endif
	};

#endif
