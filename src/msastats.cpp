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
	Aln.FromFASTAFile_PreserveCase(MSAFileName);
	f.Close();

	const uint SeqCount = Aln.GetSeqCount();
	const uint ColCount = Aln.GetColCount();
	ProgressPrefix(false);

	ProgressLog("%10u  Sequences\n", SeqCount);
	ProgressLog("%10u  Columns\n", ColCount);

	vector<uint> Ls;
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		uint L = Aln.GetUngappedSeqLength(SeqIndex);
		Ls.push_back(L);
		}
	Quarts Q;
	GetQuarts(Ls, Q);

	ProgressLog("%10.1f  Mean seq length", Q.Avg);
	ProgressLog("  min %u, median %u, max %u\n", Q.Min, Q.Med, Q.Max);

	vector<uint> GapPcts;
	uint Gap0 = 0;
	uint GapOk = 0;
	uint GapMaxPct = 0;
	float MaxGapPct = 50.0;
	uint LowerColCount = 0;
	uint UpperColCount = 0;
	uint MixedColCount = 0;
	uint DotCount = 0;
	uint DashCount = 0;
	if (optset_max_gap_fract)
		MaxGapPct = (float) (opt(max_gap_fract)*100);
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		uint NU, NL, NG, NDots, NDashes;
		Aln.GetUpperLowerGapCount(Col, NU, NL, NG, NDots, NDashes);
		DotCount += NDots;
		DashCount += NDashes;
		if (NU > 0 && NL == 0)
			++UpperColCount;
		else if (NU == 0 && NL > 0)
			++LowerColCount;
		else if (NU > 0 && NL > 0)
			++MixedColCount;
		uint GapPct = (100*NG)/SeqCount;
		if (NG == 0)
			++Gap0;
		if (GapPct <= MaxGapPct)
			++GapOk;
		GapPcts.push_back(GapPct);
		}
	GetQuarts(GapPcts, Q);

	double DashPct = GetPct(DashCount, DashCount + DotCount);

	ProgressLog("%10.1f  Mean col gap pct,", Q.Avg);
	ProgressLog(" min %u, median %u, max %u\n",
	  Q.Min, Q.Med, Q.Max);
	ProgressLog("%10u  Cols with no gaps (%.1f%% of cols)\n",
	  Gap0, GetPct(Gap0, ColCount));
	ProgressLog("%10u  Cols with <%.1f%% gaps (%.1f%% of cols)\n",
	  GapOk, MaxGapPct, GetPct(GapOk, ColCount));
	ProgressLog("%10u  Upper-case (%u lower, %u mixed)\n",
	  UpperColCount, LowerColCount, MixedColCount);
	ProgressLog("%9.1f%%  Dash gaps\n", DashPct);

	ProgressLog("\n");
	ProgressPrefix(true);
	}
