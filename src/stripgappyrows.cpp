#include "myutils.h"
#include "muscle.h"
#include "msa.h"
#include "textfile.h"
#include "seq.h"

void cmd_strip_gappy_rows()
	{
	const string &MSAFileName = opt(strip_gappy_rows);
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

	const string &OutputFileName = opt(output);
	uint DiscardCount = 0;
	FILE *f = CreateStdioFile(OutputFileName);
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		ProgressStep(SeqIndex, SeqCount, "Writing %s",
		  OutputFileName.c_str());

		const char *Label = Aln.GetSeqName(SeqIndex);
		const char *Seq = Aln.m_szSeqs[SeqIndex];
		uint GapCount = 0;
		for (uint Col = 0; Col < ColCount; ++Col)
			if (isgap(Seq[Col]))
				++GapCount;
		double GapFract = double(GapCount)/ColCount;
		if (GapFract > MaxGapFract)
			{
			++DiscardCount;
			continue;
			}
		SeqToFasta(f, (const byte *) Seq, ColCount, Label);
		}
	ProgressLog("Discarded %u / %u seqs (%.1f%%)\n",
	  DiscardCount, SeqCount, GetPct(DiscardCount, SeqCount));
	CloseStdioFile(f);
	}
