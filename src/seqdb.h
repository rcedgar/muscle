#ifndef seqdb_h
#define seqdb_h

#include "myutils.h"
#include <map>

// Minimal SeqDB host adapter for shared flat helpers (FASTA load + accessors).
class SeqDB
	{
public:
	bool m_IsAligned;
	unsigned m_ColCount;
	vector<string> m_Labels;
	vector<string> m_Seqs;
	map<string, uint> m_LabelToIndex;

public:
	SeqDB()
		{
		m_IsAligned = false;
		m_ColCount = UINT_MAX;
		}

	void Clear()
		{
		m_Seqs.clear();
		m_Labels.clear();
		m_LabelToIndex.clear();
		m_IsAligned = false;
		}

	void SetLabelToIndex();
	uint GetSeqIndex(const string &Label, bool FailOnError = true) const;
	unsigned AddSeq(const string &Label, const string &Seq);
	const string &GetSeq(unsigned SeqIndex) const;
	const byte *GetByteSeq(unsigned SeqIndex) const;
	const string &GetLabel(unsigned SeqIndex) const;
	unsigned GetSeqLength(unsigned SeqIndex) const;
	unsigned GetSeqCount() const { return SIZE(m_Seqs); }
	void FromFasta(const string &FileName, bool AllowGaps = false);
	};

#endif // seqdb_h
