#include "myutils.h"
#include "muscle.h"
#include "msa.h"
#include "textfile.h"
#include "quarts.h"

void cmd_msastats()
	{
	const string &MSAFileName = opt(msastats);
	MSA Aln;
	TextFile f(MSAFileName.c_str());
	Aln.FromFASTAFile(f);
	f.Close();

	const uint SeqCount = Aln.GetSeqCount();
	const uint ColCount = Aln.GetColCount();
	ProgressPrefix(false);

	ProgressLog("%10u  Sequences\n", SeqCount);
	ProgressLog("%10u  Columns\n", ColCount);

	vector<uint> Ls;
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		uint L = Aln.GetSeqLength(SeqIndex);
		Ls.push_back(L);
		}
	Quarts Q;
	GetQuarts(Ls, Q);

	ProgressLog("%10.1f  Mean seq length", Q.Avg);
	ProgressLog("  min %u, median %u, max %u\n", Q.Min, Q.Med, Q.Max);

	vector<uint> GapPcts;
	uint Gap0 = 0;
	uint Gap50 = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		uint GapCount = Aln.GetGapCount(Col);
		uint GapPct = (100*GapCount)/SeqCount;
		if (GapCount == 0)
			++Gap0;
		if (GapPct < 50)
			++Gap50;
		GapPcts.push_back(GapPct);
		}
	GetQuarts(GapPcts, Q);

	ProgressLog("%10.1f  Mean col gap pct,", Q.Avg);
	ProgressLog(" min %u, median %u, max %u\n",
	  Q.Min, Q.Med, Q.Max);
	ProgressLog("%10u  Cols with no gaps (%.1f%% of cols)\n",
	  Gap0, GetPct(Gap0, ColCount));
	ProgressLog("%10u  Cols with <50%% gaps (%.1f%% of cols)\n",
	  Gap50, GetPct(Gap50, ColCount));

	ProgressLog("\n");
	ProgressPrefix(true);
	}
