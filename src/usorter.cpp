#include "muscle.h"
#include "usorter.h"
#include "sort.h"

void USorter::Init()
	{
	asserta(g_AlphaSize > 0);
	if (g_Alpha == ALPHA_Amino)
		{
		m_WordLength = 3;
		m_DictSize = myipow(20, m_WordLength);
		}
	else if (g_Alpha == ALPHA_Nucleo)
		{
		m_WordLength = 8;
		m_DictSize = myipow(4, m_WordLength);
		}
	else
		asserta(false);

	m_Rows.clear();
	m_Rows.resize(m_DictSize);
	}

uint USorter::CharsToWord(const byte *Chars)
	{
	if (g_Alpha == ALPHA_Amino)
		return CharsToWord_Amino(Chars);
	else if (g_Alpha == ALPHA_Nucleo)
		return CharsToWord_Nucleo(Chars);
	asserta(false);
	return UINT_MAX;
	}

uint USorter::CharsToWord_Amino(const byte *Chars)
	{
	uint Word = 0;
	for (uint i = 0; i < m_WordLength; ++i)
		{
		char c = Chars[i];
		uint Letter = g_CharToLetter[c];
		if (Letter >= 20)
			return UINT_MAX;
		Word = Word*20 + Letter;
		}
	return Word;
	}

uint USorter::CharsToWord_Nucleo(const byte *Chars)
	{
	uint Word = 0;
	for (uint i = 0; i < m_WordLength; ++i)
		{
		char c = Chars[i];
		uint Letter = g_CharToLetter[c];
		if (Letter >= 4)
			return UINT_MAX;
		Word = Word*4 + Letter;
		}
	return Word;
	}

void USorter::AddSeq(const byte *Seq, uint L, uint SeqIndex)
	{
	asserta(g_AlphaSize > 0);
	uint Index = SIZE(m_IndexSeqIndexes);
	if (L < m_WordLength)
		return;
	const uint WordCount = L + 1 - m_WordLength;
	for (uint i = 0; i < WordCount; ++i)
		{
		uint Word = CharsToWord(Seq + i);
		if (Word < m_DictSize)
			m_Rows[Word].push_back(Index);
		}
	m_IndexSeqIndexes.push_back(SeqIndex);
	}

void USorter::SearchSeq(const byte *Seq, uint L, vector<uint> &TopSeqIndexes,
  vector<uint> &TopWordCounts)
	{
	TopSeqIndexes.clear();
	TopWordCounts.clear();

	uint IndexSize = SIZE(m_IndexSeqIndexes);
	if (IndexSize == 0)
		return;

	if (L < m_WordLength)
		return;
	const uint WordCount = L + 1 - m_WordLength;
	vector<uint> U(IndexSize, 0);
	for (uint i = 0; i < WordCount; ++i)
		{
		uint Word = CharsToWord(Seq + i);
		if (Word >= m_DictSize)
			continue;
		const vector<uint> &Row = m_Rows[Word];
		const uint n = SIZE(Row);
		for (uint j = 0; j < n; ++j)
			{
			uint TargetSeqIndex = Row[j];
			asserta(TargetSeqIndex < IndexSize);
			U[TargetSeqIndex] += 1;
			}
		}
	uint *Order = myalloc(uint, IndexSize);
	QuickSortOrderDesc(U.data(), IndexSize, Order);
	uint TopSeqIndex = Order[0];
	uint TopWordCount = U[TopSeqIndex];
	uint MinU = TopWordCount/2 - 1;
	if (MinU == 0)
		MinU = 1;
	uint LastWordCount = TopWordCount;
	for (uint i = 0; i < IndexSize; ++i)
		{
		uint Index = Order[i];
		uint WordCount = U[Index];
		asserta(WordCount <= LastWordCount);
		LastWordCount = WordCount;
		if (WordCount < MinU)
			break;
		uint SeqIndex = m_IndexSeqIndexes[Index];
		TopSeqIndexes.push_back(SeqIndex);
		TopWordCounts.push_back(WordCount);
		}
	myfree(Order);
	}

void cmd_usorter()
	{
	const string QueryFileName = opt(usorter);
	const string DBFileName = opt(db);

	MultiSequence Query;
	Query.FromFASTA(QueryFileName);

	MultiSequence DB;
	DB.FromFASTA(DBFileName);

	SetAlpha(ALPHA_Amino);

	USorter US;
	US.Init();
	const uint DBSeqCount = DB.GetSeqCount();
	for (uint DBSeqIndex = 0; DBSeqIndex < DBSeqCount; ++DBSeqIndex)
		{
		const Sequence *seq = DB.GetSequence(DBSeqIndex);
		const byte *SeqChars = (const byte *) seq->GetBytePtr();
		uint L = (uint) seq->GetLength();
		US.AddSeq(SeqChars, L, DBSeqIndex);
		}

	const uint QuerySeqCount = Query.GetSeqCount();
	for (uint QuerySeqIndex = 0; QuerySeqIndex < QuerySeqCount; ++QuerySeqIndex)
		{
		const Sequence *seq = Query.GetSequence(QuerySeqIndex);
		const byte *SeqChars = (const byte *) seq->GetBytePtr();
		uint L = (uint) seq->GetLength();

		vector<uint> TopSeqIndexes;
		vector<uint> TopWordCounts;
		US.SearchSeq(SeqChars, L, TopSeqIndexes, TopWordCounts);

		const uint n = SIZE(TopSeqIndexes);
		asserta(SIZE(TopWordCounts) == n);
		Log("\n");
		Log("Q>%s, %u hits\n", Query.GetLabel(QuerySeqIndex));
		for (uint i = 0; i < n; ++i)
			{
			uint DBSeqIndex = TopSeqIndexes[i];
			uint Count = TopWordCounts[i];
			Log("  [%4u]  %s\n", Count, DB.GetLabel(DBSeqIndex));
			}
		}
	}
