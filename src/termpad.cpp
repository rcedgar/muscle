#include "muscle.h"

void TermPad(MultiSequence &Seqs)
	{
	MultiSequence &PaddedSeqs = *new MultiSequence;
	const uint SeqCount = Seqs.GetSeqCount();
	if (SeqCount == 0)
		Die("TermPad(): empty input");
	for (uint i = 0; i < SeqCount; ++i)
		{
		Sequence *PaddedSeq = Seqs.m_Seqs[i]->TermPad();
		PaddedSeqs.AddSequence(PaddedSeq, true);
		}
	Seqs.Clear();
	Seqs = PaddedSeqs;
	}

static void GetColRange(const Sequence &seq, uint ColCount,
  uint &FirstCol, uint &LastCol)
	{
	FirstCol = UINT_MAX;
	LastCol = UINT_MAX;
	asserta(seq.GetLength() == ColCount);
	uint LeftCount = 0;
	uint RightCount = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = seq.m_CharVec[Col];
		if (c == LEFT_TERM_PAD_CHAR)
			{
			++LeftCount;
			if (LeftCount == TERM_PAD_LENGTH)
				FirstCol = Col+1;
			}
		else if (c == RIGHT_TERM_PAD_CHAR)
			{
			if (RightCount == 0)
				LastCol = Col-1;
			++RightCount;
			}
		}
	}

void DeleteTermPad(MultiSequence &Seqs)
	{
	const uint SeqCount = Seqs.GetSeqCount();
	if (SeqCount == 0)
		Die("DeleteTermPad(): empty input");
	uint ColCount = Seqs.GetColCount();

	vector<uint> FirstCols;
	vector<uint> LastCols;
	uint MinFirstCol = UINT_MAX;
	uint MaxLastCol = UINT_MAX;
	for (uint i = 0; i < SeqCount; ++i)
		{
		const Sequence &seq = *Seqs.m_Seqs[i];
		uint FirstCol, LastCol;
		GetColRange(seq, ColCount, FirstCol, LastCol);
		if (i == 0)
			{
			MinFirstCol = FirstCol;
			MaxLastCol = LastCol;
			}
		else
			{
			MinFirstCol = min(MinFirstCol, FirstCol);
			MaxLastCol = max(MaxLastCol, LastCol);
			}
		FirstCols.push_back(FirstCol);
		LastCols.push_back(LastCol);
		}

	MultiSequence &FixedSeqs = *new MultiSequence;
	for (uint i = 0; i < SeqCount; ++i)
		{
		const Sequence &seq = *Seqs.m_Seqs[i];
		Sequence &fixed_seq = *NewSequence();
		uint FirstCol = FirstCols[i];
		uint LastCol = LastCols[i];

		string str;
		for (uint Col = MinFirstCol; Col < FirstCol; ++Col)
			str += '-';
		for (uint Col = FirstCol; Col <= LastCol; ++Col)
			str += seq.m_CharVec[Col];
		for (uint Col = LastCol+1; Col <= MaxLastCol; ++Col)
			str += '-';
		fixed_seq.FromString(seq.m_Label, str);
		FixedSeqs.AddSequence(&fixed_seq, true);
		}

	Seqs.Clear();
	Seqs = FixedSeqs;
	}
