#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <cctype>
#include <cstdlib>
#include "FileBuffer.h"

// data uses one-based indexing, data[0] set to '@'.
class Sequence
	{
private:
// Global input MulitSquence index
	uint m_GSI = UINT_MAX;

// Sparse matrix index
	uint m_SMI = UINT_MAX;

public:
	string label;
	vector<char>* data;

public:
	Sequence()
		{
		label = "~";
		data = 0;
		m_GSI = UINT_MAX;
		m_SMI = UINT_MAX;
		}

	bool FromFileBuffer(FileBuffer& infile, bool stripGaps = false);

	void Create(vector<char>* data, string label, uint GSI, uint SMI);

	~Sequence()
		{
		if (data)
			{
			delete data;
			data = NULL;
			}
		}

	void FromString(const string &Label, const string &Seq);
	void Copy(const Sequence &rhs);

	void InitData()
		{
		asserta(data == 0);
		data = new vector<char>;
		asserta(data);
		(*(data)).push_back('@');
		}

	void AppendChar(char c)
		{
		(*(data)).push_back(c);
		}

	const string &GetLabel() const
		{
		return label;
		}

	const char *GetLabelCStr() const
		{
		return label.c_str();
		}

// One-based
	char *GetCharPtr1()
		{
		assert(data);
		return data->data();
		}

	const byte *GetBytePtr() const
		{
		assert(data);
		const char *CharPtr = data->data() + 1;
		const byte *BytePtr = (const byte *) CharPtr;
		return BytePtr;
		}

// Chars stored with one-based indexing.
	char GetPosition(int i) const
		{
		assert(data != 0);
		assert(i >= 1 && i < data->size());
		return (*data)[i];
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
		label = Label;
		}

	void SetSMI(uint SMI)
		{
		m_SMI = SMI;
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
		if (data == 0)
			return 0;
		uint n = (uint) data->size();
		asserta(n > 0);
		uint L = n - 1;
		return L;
		}

	void WriteMFA(FILE *f) const;
	Sequence* Clone() const;
	Sequence* AddGaps(const vector<char>* alignment, char id) const;
	Sequence* DeleteGaps() const;
	void GetPosToCol_OneBased(vector<uint> &PosToCol) const;
	void GetPosToCol(vector<uint> &PosToCol) const;
	void GetColToPos(vector<uint> &ColToPos) const;
	void LogMe() const;
	};

#endif
