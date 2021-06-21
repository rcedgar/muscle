#ifndef MULTISEQUENCE_H
#define MULTISEQUENCE_H

#include <cctype>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <set>
#include "sequence.h"

class MultiSequence
	{
public:
	vector<Sequence*>* sequences = 0;

public:
	MultiSequence()
		{
		sequences = 0;
		}

	MultiSequence(FileBuffer& infile)
		{
		sequences = 0;
		LoadMFA(infile);
		}

	MultiSequence(const string& filename)
		{
		sequences = 0;
		LoadMFA(filename);
		}

	~MultiSequence();
	void FreeNonOwner();
	void Copy(const MultiSequence &rhs);
	void LoadMFA(const string& filename, bool stripGaps = false);
	void LoadMFA(FileBuffer& infile, bool stripGaps = false);
	void FromFASTA(const string& filename, bool stripGaps = false)
		{
		LoadMFA(filename, stripGaps);
		}
	uint GetSeqIndex(const string &Label, bool FailOnError = true) const;

	void AddSequence(Sequence* sequence)
		{
		assert(sequence);

		// add sequence
		if (!sequences) sequences = new vector<Sequence*>;
		sequences->push_back(sequence);
		}

	void WriteMFA(FILE *f) const
		{
		if (f == 0)
			return;
		if (!sequences)
			return;

		for (vector<Sequence*>::iterator iter = sequences->begin();
		  iter != sequences->end(); ++iter)
			{
			const Sequence *Seq = *iter;
			asserta(Seq != 0);
			Seq->WriteMFA(f);
			}
		}

	void WriteMFA(const string &FileName) const
		{
		if (FileName.empty())
			return;
		FILE *f = CreateStdioFile(FileName);
		WriteMFA(f);
		CloseStdioFile(f);
		}

	Sequence* GetSequence(int i)
		{
		assert(sequences);
		assert(0 <= i && i < (int)sequences->size());
		return (*sequences)[i];
		}

	const Sequence* GetSequence(int i) const
		{
		assert(sequences);
		assert(0 <= i && i < (int)sequences->size());
		return (*sequences)[i];
		}

	int GetNumSequences() const
		{
		if (!sequences)
			return 0;
		return (int)sequences->size();
		}

	uint GetSeqCount() const
		{
		if (sequences == 0)
			return 0;
		return SIZE(*sequences);
		}

	uint GetChar(uint SeqIndex, uint ZeroBasedPos) const
		{
		const Sequence *seq = GetSequence((int) SeqIndex);
		char c = seq->GetPosition(int(ZeroBasedPos)+1);
		return c;
		}

	MultiSequence* Project(const set<int>& indices);

	bool IsAligned() const;
	uint GetColCount() const;

	const string &GetLabelStr(uint SeqIndex) const
		{
		const Sequence *seq = GetSequence((int) SeqIndex);
		return seq->label;
		}

	const char *GetLabel(uint SeqIndex) const
		{
		const Sequence *seq = GetSequence((int) SeqIndex);
		return seq->label.c_str();
		}

	const byte *GetByteSeq(uint SeqIndex, uint &L) const
		{
		const Sequence *seq = GetSequence((int) SeqIndex);
		const byte *ByteSeq = (const byte *) seq->data->data() + 1;
		L = (uint) seq->GetLength();
		return ByteSeq;
		}
	bool GuessIsNucleo() const;
	void LogGSIs(const char *Msg = 0) const;
	void GetLengthOrder(vector<uint> &SeqIndexes) const;
	uint GetLength(uint SeqIndex) const;
	uint GetGSI(uint SeqIndex) const;
	};

#endif
