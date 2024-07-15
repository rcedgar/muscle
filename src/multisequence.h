#pragma once

#include <set>
#include "sequence.h"


class MSA;
class Mega;

class MultiSequence
	{
public:
	vector<const Sequence *> m_Seqs;
	vector<bool> m_Owners;
	bool m_DupeLabelsOk = false;

public:
	static uint m_NewCount;
	static uint m_DeleteCount;

public:
	MultiSequence()
		{
		//Clear();
#pragma omp critical
		++m_NewCount;
		}

	MultiSequence(FileBuffer& infile)
		{
		LoadMFA(infile);
#pragma omp critical
		++m_NewCount;
		}

	MultiSequence(const string& filename)
		{
		LoadMFA(filename);
#pragma omp critical
		++m_NewCount;
		}

	void Clear();
	~MultiSequence()
		{ 
		Clear();
#pragma omp critical
		++m_DeleteCount;
		}

	void FromStrings(const vector<string> &Labels, const vector<string> &Seqs);
	void Copy(const MultiSequence &rhs);
        
    void LoadMega(Mega & MM);
	void LoadMFA(const string& filename, bool stripGaps = false);
	void LoadMFA(FileBuffer& infile, bool stripGaps = false);
	void FromFASTA(const string& filename, bool stripGaps = false)
		{
		LoadMFA(filename, stripGaps);
		}
	uint GetSeqIndex(const string &Label, bool FailOnError = true) const;
	void ToMSA(MSA &msa) const;

	void AddSequence(const Sequence *sequence, bool Owner)
		{
		m_Seqs.push_back(sequence);
		m_Owners.push_back(Owner);
		}

	void ToFasta(FILE *f) const
		{
		WriteMFA(f);
		}

	void ToFasta(const string &FileName) const
		{
		WriteMFA(FileName);
		}

	void WriteMFA(FILE *f) const
		{
		if (f == 0)
			return;

		for (uint i = 0; i < SIZE(m_Seqs); ++i)
			{
			const Sequence *Seq = m_Seqs[i];
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

	const Sequence *GetSequence(int i) const
		{
		return m_Seqs[i];
		}

	int GetNumSequences() const
		{
		return (int) SIZE(m_Seqs);
		}

	uint GetSeqCount() const
		{
		return SIZE(m_Seqs);
		}

	double GetMeanSeqLength() const;
	uint GetMaxSeqLength() const;
	uint GetMinSeqLength() const;

	uint GetChar(uint SeqIndex, uint ZeroBasedPos) const
		{
		const Sequence *seq = GetSequence((int) SeqIndex);
		char c = seq->GetChar(int(ZeroBasedPos));
		return c;
		}

	MultiSequence* Project(const set<int>& indices);
	MultiSequence* Project(const vector<uint> &SeqIndexes);

	bool IsAligned() const;
	uint GetColCount() const;

	const string &GetLabelStr(uint SeqIndex) const
		{
		const Sequence *seq = GetSequence((int) SeqIndex);
		return seq->m_Label;
		}

	void GetSeqStr(uint SeqIndex, string &Seq) const
		{
		const Sequence *seq = GetSequence((int) SeqIndex);
		seq->GetSeqAsString(Seq);
		}

	const char *GetLabel(uint SeqIndex) const
		{
		const Sequence *seq = GetSequence((int) SeqIndex);
		return seq->m_Label.c_str();
		}

	const byte *GetByteSeq(uint SeqIndex, uint &L) const
		{
		const Sequence *seq = GetSequence((int) SeqIndex);
		const byte *ByteSeq = (const byte *) seq->m_CharVec.data();
		L = (uint) seq->GetLength();
		return ByteSeq;
		}

	const byte *GetBytePtr(uint SeqIndex) const
		{
		const Sequence *seq = GetSequence((int) SeqIndex);
		const byte *BytePtr = seq->GetBytePtr();
		return BytePtr;
		}

	const char *GetCharPtr(uint SeqIndex) const
		{
		const Sequence *seq = GetSequence((int) SeqIndex);
		const char *CharPtr = seq->GetCharPtr();
		return CharPtr;
		}

	bool GuessIsNucleo() const;
	void GetLengthOrder(vector<uint> &SeqIndexes) const;
	uint GetSeqLength(uint SeqIndex) const;
	void LogMe() const;
	bool ColIsAllGaps(uint Col) const;
	};
