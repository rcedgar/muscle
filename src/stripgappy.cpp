#include "myutils.h"
#include "muscle.h"
#include "msa.h"
#include "textfile.h"
#include "seq.h"

void cmd_strip_gappy()
	{
	const string &MSAFileName = opt(strip_gappy);
	double MaxGapFract = 0.5;
	double MaxGapFractRow = 0.5;
	if (optset_max_gap_fract)
		MaxGapFract = opt(max_gap_fract);
	if (optset_max_gap_fract_row)
		MaxGapFractRow = opt(max_gap_fract_row);

	MSA Aln;
	TextFile TF(MSAFileName.c_str());
	Aln.FromFASTAFile(TF);
	TF.Close();

	const uint SeqCount = Aln.GetSeqCount();
	const uint ColCount = Aln.GetColCount();
	ProgressLog("%u seqs, %u cols, max gaps %.4f, %.3g\n",
	  SeqCount, ColCount, MaxGapFract, MaxGapFractRow);

	vector<uint> KeepCols;
	uint DiscardColCount = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		uint GapCount = Aln.GetGapCount(Col);
		double GapFract = double(GapCount)/double(SeqCount);
		if (GapFract <= MaxGapFract)
			KeepCols.push_back(Col);
		else
			++DiscardColCount;
		}

	uint NewColCount = SIZE(KeepCols);
	asserta(NewColCount > 0);
	//ProgressLog("Keeping %u cols (%.1f%%)\n",
	//  NewColCount, GetPct(NewColCount, ColCount));

	const string &OutputFileName = opt(output);
	FILE *f = CreateStdioFile(OutputFileName);
	byte *NewSeq = myalloc(byte, NewColCount);
	uint DiscardRowCount = 0;
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		ProgressStep(SeqIndex, SeqCount, "Writing %s",
		  OutputFileName.c_str());
		const char *Label = Aln.GetSeqName(SeqIndex);
		uint RowGapCount = 0;
		for (uint i = 0; i < NewColCount; ++i)
			{
			uint Col = KeepCols[i];
			char c = Aln.GetChar(SeqIndex, Col);
			NewSeq[i] = c;
			if (c == '-')
				++RowGapCount;
			}
		double GapFractRow = double(RowGapCount)/NewColCount;
		if (GapFractRow > MaxGapFractRow)
			{
			++DiscardRowCount;
			continue;
			}
		SeqToFasta(f, NewSeq, NewColCount, Label);
		}
	ProgressLog("%u cols, %u rows discarded (max gaps %.4g, %.3g)\n",
	  DiscardColCount, DiscardRowCount, MaxGapFract, MaxGapFractRow);
	CloseStdioFile(f);
	}
