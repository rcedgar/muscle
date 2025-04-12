#include "muscle.h"
#include "msa.h"
#include "alpha.h"

const double MAX_GAP_FRACT = 0.1;

void GetInsertLoHis(const vector<bool> &IsAlignedVec,
  vector<uint> &Los, vector<uint> &His);

static bool GetMSAColIsAligned(const MultiSequence &Aln, uint Col,
							 double MaxGapFract)
	{
	const uint SeqCount = Aln.GetSeqCount();
	if (SeqCount == 0)
		return false;

	uint GapCount = 0;
	uint UpperCount = 0;
	uint LowerCount = 0;
	for (uint i = 0; i < SeqCount; ++i)
		{
		char c = Aln.GetChar(i, Col);
		if (isgap(c))
			++GapCount;
		}
	double GapFract = double(GapCount)/(SeqCount+0.01);
	return GapFract < MaxGapFract;
	}

static uint GetMSAColAlignedVec(const MultiSequence &Aln, vector<bool> &GappyVec)
	{
	GappyVec.clear();
	const uint ColCount = Aln.GetColCount();
	uint n = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		bool gappy = GetMSAColIsAligned(Aln, Col, MAX_GAP_FRACT);
		if (gappy)
			++n;
		GappyVec.push_back(gappy);
		}
	return n;
	}

MultiSequence *SqueezeGappyCols(const MultiSequence &Aln)
	{
	const uint SeqCount = Aln.GetSeqCount();
	const uint ColCount = Aln.GetColCount();

	vector<string> UngappedSeqs;
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		string UngappedSeq;
		Aln.GetUngappedSeqStr(SeqIndex, UngappedSeq);
		UngappedSeqs.push_back(UngappedSeq);
		}

	vector<bool> IsAlignedVec;
	uint n = GetMSAColAlignedVec(Aln, IsAlignedVec);
	ProgressLog("%u / %u gappy cols\n", ColCount - n, ColCount);

	vector<uint> Los;
	vector<uint> His;
	GetInsertLoHis(IsAlignedVec, Los, His);
	const uint RangeCount = SIZE(Los);

	if (RangeCount == 0)
		{
		MultiSequence *S = new MultiSequence;
		S->Copy(Aln);
		return S;
		}

	asserta(SIZE(His) == RangeCount);
	vector<string> NewRows(SeqCount);
	uint PrevHi = UINT_MAX;
	for (uint RangeIndex = 0; RangeIndex < RangeCount; ++RangeIndex)
		{
		uint Lo = Los[RangeIndex];
		uint Hi = His[RangeIndex];
		asserta(Lo <= Hi);
		if (RangeIndex > 0)
			asserta(Lo > PrevHi);

		uint From = (PrevHi == UINT_MAX ? 0 : PrevHi + 1);
		for (uint Col = From; Col < Lo; ++Col)
			for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
				{
				char c =  Aln.GetChar(SeqIndex, Col);
				NewRows[SeqIndex] += c;
				}

		uint MaxInsL = 0;
		vector<string> Inserts;
		for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
			{
			string Insert;
			for (uint Col = Lo; Col <= Hi; ++Col)
				{
				char c = Aln.GetChar(SeqIndex, Col);
				if (!isgap(c))
					Insert += c;
				}
			uint L = SIZE(Insert);
			if (L > MaxInsL)
				MaxInsL = L;
			Inserts.push_back(Insert);
			}
		for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
			{
			const string &Insert = Inserts[SeqIndex];
			uint InsL = SIZE(Insert);
			uint n = 0;
			uint Dots = MaxInsL - InsL;
			uint Dots1 = Dots/2;
			if (From == 0)
				Dots1 = Dots;
			else if (RangeIndex + 1 == RangeCount)
				Dots1 = 0;
			while (n < Dots1)
				{
				NewRows[SeqIndex] += "-";
				++n;
				}
			NewRows[SeqIndex] += Insert;
			asserta(InsL <= MaxInsL);
			while (n < Dots)
				{
				NewRows[SeqIndex] += "-";
				++n;
				}
			}
		for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
			asserta(SIZE(NewRows[SeqIndex]) == SIZE(NewRows[0]));

		PrevHi = Hi;
		}

	uint From = (PrevHi == UINT_MAX ? 0 : PrevHi + 1);
	for (uint Col = From; Col < ColCount; ++Col)
		for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
			NewRows[SeqIndex] += Aln.GetChar(SeqIndex, Col);
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		asserta(SIZE(NewRows[SeqIndex]) == SIZE(NewRows[0]));

	vector<string> Labels;
	const uint NewColCount = SIZE(NewRows[0]);
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		string Label = (string) Aln.GetLabel(SeqIndex);
		Labels.push_back(Label);
		}

	MultiSequence *S = new MultiSequence;
	S->FromStrings2(Labels, NewRows);

	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		string UngappedSeq;
		S->GetUngappedSeqStr(SeqIndex, UngappedSeq);
		if (UngappedSeq != UngappedSeqs[SeqIndex])
			{
			Log("\n");
			Log(">%s\n", Aln.GetLabel(SeqIndex));
			Log("%s\n", UngappedSeq.c_str());
			Log("%s\n", UngappedSeqs[SeqIndex].c_str());
			Die("UngappedSeq != UngappedSeqs[SeqIndex]");
			}
		}
	ProgressLog("%u squeezed cols\n", S->GetColCount());
	return S;
	}

void cmd_squeeze_gappy()
	{
	const string &InputFileName = g_Arg1;

	MultiSequence InputMSA;
	InputMSA.FromFASTA(InputFileName);

	MultiSequence *OutMSA = SqueezeGappyCols(InputMSA);

	OutMSA->ToFasta(opt(output));
	}
