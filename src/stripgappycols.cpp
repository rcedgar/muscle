#include "myutils.h"
#include "muscle.h"
#include "msa.h"
#include "textfile.h"
#include "seq.h"

void cmd_strip_gappy_cols()
	{
	const string &MSAFileName = opt(strip_gappy_cols);
	double MaxGapFract = 0.5;
	if (optset_max_gap_fract)
		MaxGapFract = opt(max_gap_fract);

	MSA Aln;
	TextFile TF(MSAFileName.c_str());
	Aln.FromFASTAFile(TF);
	TF.Close();

	const uint SeqCount = Aln.GetSeqCount();
	const uint ColCount = Aln.GetColCount();
	ProgressLog("%u seqs, %u cols, max gaps %.4f\n",
	  SeqCount, ColCount, MaxGapFract);

	vector<uint> KeepCols;
	uint GappyCount = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		uint GapCount = Aln.GetGapCount(Col);
		double GapFract = double(GapCount)/double(SeqCount);
		if (GapFract <= MaxGapFract)
			KeepCols.push_back(Col);
		else
			++GappyCount;
		}

	uint NewColCount = SIZE(KeepCols);
	ProgressLog("Keeping %u cols (%.1f%%)\n",
	  NewColCount, GetPct(NewColCount, ColCount));
	asserta(NewColCount > 0);

	const string &OutputFileName = opt(output);
	FILE *f = CreateStdioFile(OutputFileName);
	byte *NewSeq = myalloc(byte, NewColCount);
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		ProgressStep(SeqIndex, SeqCount, "Writing %s",
		  OutputFileName.c_str());
		const char *Label = Aln.GetSeqName(SeqIndex);
		for (uint i = 0; i < NewColCount; ++i)
			{
			uint Col = KeepCols[i];
			char c = Aln.GetChar(SeqIndex, Col);
			NewSeq[i] = c;
			}
		SeqToFasta(f, NewSeq, NewColCount, Label);
		}
	CloseStdioFile(f);
	}
