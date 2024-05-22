#include "muscle.h"
#include "msa.h"
#include "alpha.h"

MSA *SqueezeInserts(const MSA &Aln, ptr_GetMSAColIsAligned pFn);
bool GetMSAColIsAligned(const MSA &Aln, uint Col);

const uint MINCORECOLS = 10;

void cmd_mustang_core()
	{
	const string &InputFileName = g_Arg1;

	MSA InputMSA;
	InputMSA.FromFASTAFile_PreserveCase(InputFileName);
	const uint SeqCount = InputMSA.GetSeqCount();
	const uint ColCount = InputMSA.GetColCount();

	MSA *TmpMSA = new MSA;
	TmpMSA->Copy(InputMSA);

	uint AlignedColCount = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		//bool Aligned = GetMSAColIsAligned(InputMSA, Col);
		bool Aligned = !InputMSA.ColumnHasGap(Col);
		if (Aligned)
			++AlignedColCount;
		for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
			{
			char c = InputMSA.GetChar(SeqIndex, Col);
			if (Aligned)
				{
				if (isgap(c))
					c = '-';
				else if (isalpha(c))
					c = toupper(c);
				}
			else
				{
				if (isgap(c))
					c = '.';
				else if (isalpha(c))
					c = tolower(c);
				}
			TmpMSA->SetChar(SeqIndex, Col, c);
			}
		}

	MSA *OutMSA = SqueezeInserts(*TmpMSA, GetMSAColIsAligned);
	asserta(OutMSA->GetSeqCount() == SeqCount);
	const uint NewColCount = OutMSA->GetColCount();

	if (AlignedColCount < MINCORECOLS)
		{
		Warning("%u aligned cols < %u", AlignedColCount, MINCORECOLS);
		return;
		}

	FILE *f = CreateStdioFile(opt(output));
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const char *Seq = OutMSA->GetSeqCharPtr(SeqIndex);
		string Label;
		OutMSA->GetSeqLabel(SeqIndex, Label);
		if (EndsWith(Label, ".pdb"))
			Label = Label.substr(0, Label.size() - 4);
		SeqToFasta(f, Seq, NewColCount, Label.c_str());
		}
	CloseStdioFile(f);
	}
