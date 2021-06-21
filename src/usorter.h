#pragma once

class USorter
	{
public:
	const MultiSequence m_MFA;
	vector<vector<uint> > m_Rows;
	vector<uint> m_IndexSeqIndexes;
	uint m_WordLength = 3;
	uint m_DictSize = myipow(20, 3);

public:
	void Init();
	void AddSeq(const byte *Seq, uint L, uint SeqIndex);
	uint CharsToWord(const byte *Chars);
	void SearchSeq(const byte *Seq, uint L,
	  vector<uint> &TopSeqIndexes, vector<uint> &TopWordCounts);
	};
