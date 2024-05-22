#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <cctype>
#include <cstdlib>
#include "filebuffer.h"

class Sequence
	{
public:
// Global input MulitSquence index
	uint m_GSI = UINT_MAX;

// Sparse matrix index
	uint m_SMI = UINT_MAX;

	string m_Label;
	vector<char> m_CharVec;

private:
	Sequence()
		{
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

	void InitData()
		{
		m_CharVec.clear();
		//m_CharVec.push_back('@');
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

	const char *GetCharPtr() const
		{
		const char *CharPtr = m_CharVec.data();
		return CharPtr;
		}

	const byte *GetBytePtr() const
		{
		const char *CharPtr = m_CharVec.data();
		const byte *BytePtr = (const byte *) CharPtr;
		return BytePtr;
		}

	char GetChar(uint i) const
		{
		assert(i >= 0 && i < m_CharVec.size());
		return m_CharVec[i];
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
		return SIZE(m_CharVec);
		}

	void WriteMFA(FILE *f) const;
	Sequence* Clone() const;
	Sequence* AddGapsPath(const string &Path, char id) const;
	Sequence* CopyDeleteGaps() const;
	void GetPosToCol(vector<uint> &PosToCol) const;
	void GetColToPos(vector<uint> &ColToPos) const;
	void LogMe() const;
	void GetSeqAsString(string &Seq) const;
public:
	static Sequence *_NewSequence();
#define NewSequence()	Sequence::_NewSequence()
	static void _DeleteSequence(const Sequence *Seq);
#define DeleteSequence(s)	Sequence::_DeleteSequence((s))

	static void LogNewDeleteCounts();
	};

#endif
