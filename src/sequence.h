#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <cctype>
#include <cstdlib>
#include "filebuffer.h"

//static const char LEFT_TERM_PAD_CHAR = '<';
//static const char RIGHT_TERM_PAD_CHAR = '>';
//static const uint TERM_PAD_LENGTH = 3;

class Sequence
	{
public:
	string m_Label;
	vector<char> m_CharVec;

private:
	Sequence() {}
	~Sequence() {}

public:
	bool FromFileBuffer(FileBuffer& infile, bool stripGaps = false);
	void FromString(const string &Label, const string &Seq);

	void InitData()
		{
		m_CharVec.clear();
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

	void OverwriteLabel(const string &Label)
		{
		m_Label = Label;
		}

	uint GetLength() const
		{
		return SIZE(m_CharVec);
		}

	void WriteMFA(FILE *f) const;
	Sequence *Clone() const;
	Sequence *AddGapsPath(const string &Path, char id) const;
	Sequence *CopyDeleteGaps() const;
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
