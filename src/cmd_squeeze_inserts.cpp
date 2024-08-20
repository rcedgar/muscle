#include "muscle.h"
#include "msa.h"
#include "alpha.h"

void GetInsertLoHis(const vector<bool> &Als,
  vector<uint> &Los, vector<uint> &His)
	{
	Los.clear();
	His.clear();
	uint Start = UINT_MAX;
	const uint ColCount = SIZE(Als);
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		if (Als[Col])
			{
			if (Start != UINT_MAX)
				{
				Los.push_back(Start);
				His.push_back(Col-1);
				Start = UINT_MAX;
				}
			}
		else
			{
			if (Start == UINT_MAX)
				Start = Col;
			}
		}
	if (Start != UINT_MAX)
		{
		Los.push_back(Start);
		His.push_back(ColCount-1);
		}
	}

bool GetMSAColIsAligned(const MSA &Aln, uint Col)
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
		else if (isupper(c))
			++UpperCount;
		else if (islower(c))
			++LowerCount;
		else
			{
			Aln.LogMe();
			Die("Unexpected sequence char '%c'", c);
			}
		}
	if (optset_max_gap_fract)
		{
		double GapFract = double(GapCount)/SeqCount;
		return GapFract <= opt(max_gap_fract);
		}
	else
		{
		if (UpperCount > 0 && LowerCount > 0)
			Die("Mixed-case col");
		return UpperCount > 0;
		}
	}

void GetMSAColAlignedVec(const MSA &Aln,
  vector<bool> &AlignedVec, ptr_GetMSAColIsAligned pFn)
	{
	AlignedVec.clear();
	const uint ColCount = Aln.GetColCount();
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		//bool Al = GetMSAColIsAligned(Aln, Col);
		bool Al = pFn(Aln, Col);
		AlignedVec.push_back(Al);
		}
	}

MSA *SqueezeInserts(const MSA &Aln, ptr_GetMSAColIsAligned pFn)
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

	vector<bool> Als;
	GetMSAColAlignedVec(Aln, Als, pFn);

	vector<uint> Los;
	vector<uint> His;
	GetInsertLoHis(Als, Los, His);
	const uint RangeCount = SIZE(Los);

	if (RangeCount == 0)
		{
		MSA *S = new MSA;
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
				if (Als[Col])
					c = toupper(c);
				else
					c = tolower(c);
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
					Insert += tolower(c);
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
				NewRows[SeqIndex] += ".";
				++n;
				}
			NewRows[SeqIndex] += Insert;
			asserta(InsL <= MaxInsL);
			while (n < Dots)
				{
				NewRows[SeqIndex] += ".";
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

	MSA *S = new MSA;
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
	return S;
	}

void cmd_squeeze_inserts()
	{
	const string &InputFileName = g_Arg1;

	MSA InputMSA;
	InputMSA.FromFASTAFile_PreserveCase(InputFileName);

	MSA *OutMSA = SqueezeInserts(InputMSA, GetMSAColIsAligned);

	OutMSA->ToFASTAFile(opt(output));
	}
