#include "muscle.h"
#include "ensemble.h"
#include "sort.h"

void cmd_efa_bestcols()
	{
	const string EfaFileName = opt(efa_bestcols);
	const string OutputFileName = opt(output);

	double MinConf = 1.0;
	if (optset_minconf)
		MinConf = opt(minconf);

	double MaxGapFract = 0.5;
	if (optset_max_gap_fract)
		MaxGapFract = opt(max_gap_fract);
	asserta(MaxGapFract >= 0 && MaxGapFract <= 1.0);

	uint MaxCols = UINT_MAX;
	if (optset_maxcols)
		MaxCols = opt(maxcols);

	Ensemble E;
	E.FromFile(EfaFileName);

	vector<double> Confs;
	const uint UniqueIxCount = SIZE(E.m_UniqueIxToIxs);
	vector<uint> UniqueIxs;
	for (uint UniqueIx = 0; UniqueIx < UniqueIxCount; ++UniqueIx)
		{
		uint n = SIZE(Confs);
		double Pct = GetPct(n, UniqueIxCount);
		ProgressStep(UniqueIx, UniqueIxCount,
		  "%u cols (%.1f%%) conf >= %.3g, gaps <= %.3g",
		  n, Pct, MinConf, MaxGapFract);

		double Conf = E.GetConf(UniqueIx);
		if (Conf < MinConf)
			continue;
		uint Ix = E.m_UniqueIxs[UniqueIx];
		double GapFract = E.GetGapFract(Ix);
		if (GapFract > MaxGapFract)
			continue;

		UniqueIxs.push_back(UniqueIx);
		Confs.push_back(Conf);
		}

	const uint M = SIZE(Confs);
	vector<uint> Order(M);
	QuickSortOrderDesc<double>(Confs.data(), M, Order.data());

	vector<uint> BestUniqueIxs;
	const uint N = min(SIZE(Order), MaxCols);
	for (uint i = 0; i < N; ++i)
		{
		uint UniqueIx = UniqueIxs[Order[i]];
		BestUniqueIxs.push_back(UniqueIx);
		}

	MSA RepAln;
	E.MakeResampledMSA(BestUniqueIxs, RepAln);
	RepAln.ToFASTAFile(OutputFileName);
	}
