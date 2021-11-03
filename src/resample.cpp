#include "muscle.h"
#include "ensemble.h"

void cmd_resample()
	{
	const string &FileName = opt(resample);
	const string &OutputPattern = opt(output);
	if (OutputPattern.empty())
		Die("Must set -output");

	double MaxGapFract = 0.5;
	if (optset_max_gap_fract)
		MaxGapFract = opt(max_gap_fract);

	double MinConf = 0.5;
	if (optset_minconf)
		MinConf = opt(minconf);

	uint ReplicateCount = 100;
	if (optset_replicates)
		ReplicateCount = opt(replicates);

	Ensemble E;
	E.FromFile(FileName);

	uint SiteCount = E.GetMedianHiQualColCount(MaxGapFract, MinConf);
	if (SiteCount == 0)
		Die("All columns low qual (max fract %.3g, min conf %.3g)",
		  MaxGapFract, MinConf);
	ProgressLog("Site count %u\n", SiteCount);
	if (SiteCount < 20)
		Warning("Very low hi qual site count");

	vector<uint> NonGappyUniqueIxs;
	E.GetHiQualUniqueIxs(MaxGapFract, MinConf, NonGappyUniqueIxs);
	const uint N = SIZE(NonGappyUniqueIxs);

	bool OutputWildCard = (OutputPattern.find('@') != string::npos);
	FILE *fOut = 0;
	if (!OutputWildCard)
		fOut = CreateStdioFile(OutputPattern);
	for (uint RepIndex = 0; RepIndex < ReplicateCount; ++RepIndex)
		{
		if (ReplicateCount > 1)
			ProgressStep(RepIndex, ReplicateCount, "Resampling");

		vector<uint> ResampledUniqueIxs;
		for (uint i = 0; i < SiteCount; ++i)
			{
			uint r = randu32()%N;
			uint UniqueIx = NonGappyUniqueIxs[r];
			ResampledUniqueIxs.push_back(UniqueIx);
			}

		MSA RepAln;
		E.MakeResampledMSA(ResampledUniqueIxs, RepAln);
		if (OutputWildCard)
			{
			string OutputFileName;
			MakeReplicateFileName_N(OutputPattern, RepIndex+1, OutputFileName);
			fOut = CreateStdioFile(OutputFileName);
			}
		else
			Pf(fOut, "<resampled.%u\n", RepIndex+1);
		RepAln.ToFASTAFile(fOut);
		if (OutputWildCard)
			CloseStdioFile(fOut);
		}

	if (!OutputWildCard)
		CloseStdioFile(fOut);
	}
