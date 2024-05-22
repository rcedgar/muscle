#include "muscle.h"

void cmd_make_a2m_refseq()
	{
	const string &InputFileName = opt(make_a2m_refseq);
	MSA Aln;
	Progress("Reading %s ...", InputFileName.c_str());
	Aln.FromFASTAFile(InputFileName);
	Progress("done.\n");

	FILE *fOut = CreateStdioFile(opt(output));

	const uint SeqCount = Aln.GetSeqCount();
	const uint ColCount = Aln.GetColCount();
	asserta(SeqCount > 0);
	asserta(ColCount > 0);

	uint RefSeqIndex = 0;
	string RefLabel;
	if (optset_label)
		{
		RefLabel = opt(label);
		RefSeqIndex = Aln.GetSeqIndex(RefLabel);
		}
	else
		Aln.GetSeqLabel(0, RefLabel);

	vector<uint> RefPosToCol;
	vector<uint> RefColToPos;
	Aln.GetPosToCol(RefSeqIndex, RefPosToCol);
	Aln.GetColToPos(RefSeqIndex, RefColToPos);

	const uint RL = SIZE(RefPosToCol);

	vector<bool> IsInserts;
	IsInserts.resize(ColCount, false);

	uint FirstCol = RefPosToCol[0];
	for (uint Col = 0; Col < FirstCol; ++Col)
		IsInserts[Col] = true;

	for (uint RefPos = 1; RefPos < RL; ++RefPos)
		{
		uint PrevCol = RefPosToCol[RefPos-1];
		uint ThisCol = RefPosToCol[RefPos];
		asserta(PrevCol < ThisCol);

		for (uint Col = PrevCol + 1; Col < ThisCol; ++Col)
			IsInserts[Col] = true;
		}

	uint LastCol = RefPosToCol[RL-1];
	for (uint Col = LastCol + 1; Col < ColCount; ++Col)
		IsInserts[Col] = true;

	for (uint Col = 0; Col < ColCount; ++Col)
		{
		bool IsMatch = (RefColToPos[Col] != UINT_MAX);
		bool IsInsert = IsInserts[Col];
		asserta(int(IsMatch) + int(IsInsert) == 1);
		}

	for (uint i = 0; i < SeqCount; ++i)
		{
		uint SeqIndex = i;
		if (i == 0)
			SeqIndex = RefSeqIndex;
		else if (i == RefSeqIndex)
			SeqIndex = 0;
		ProgressStep(i, SeqCount, "Converting");
		Pf(fOut, ">%s\n", Aln.GetSeqName(SeqIndex));

		string SeqStr;
		const char *SeqCharPtr = Aln.GetSeqCharPtr(SeqIndex);
		for (uint Col = 0; Col < ColCount; ++Col)
			{
			char c = SeqCharPtr[Col];

			bool IsMatch = (RefColToPos[Col] != UINT_MAX);
			bool IsInsert = IsInserts[Col];
			
			extern bool g_FASTA_AllowDigits;
			if (IsMatch)
				{
				asserta(!IsInsert);
				if (c == '.' || c == '-')
					c = '-';
				else if (isalpha(c))
					c = toupper(c);
				else if (g_FASTA_AllowDigits && isdigit(c))
					;
				else
					Die("Bad char 0x%02x", c);
				}
			else if (IsInsert)
				{
				asserta(!IsMatch);
				if (c == '.' || c == '-')
					continue;
				c = tolower(c);
				}
			else
				asserta(false);
			SeqStr += c;
			}

		Pf(fOut, "%s\n", SeqStr.c_str());
		}

	CloseStdioFile(fOut);
	}
