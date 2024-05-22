#include "muscle.h"
#include "bench.h"
#include "sweeper.h"
#include "m3alnparams.h"

void ParseGridSpec(const string &Spec,
  vector<string> &Names,
  vector<float> &Goods,
  vector<float> &Los,
  vector<float> &His,
  vector<uint> &Sizes);

static Bench *g_Bench;
static double g_TopScore = 0;
static double g_TopQ = 0;
static double g_TopTC = 0;
static float g_SubstMx_Letter[20][20];
static M3AlnParams g_AP;

static void SweeperGetScore(const Sweeper &S, const vector<float> &ParamValues,
  double &Q, double &TC)
	{
	printf("%s  ", GetProgressPrefixCStr());
	asserta(SIZE(ParamValues) == S.m_ParamCount);
	float GapOpen = FLT_MAX;
	float Center = FLT_MAX;
	for (uint ParamIndex = 0; ParamIndex < S.m_ParamCount; ++ParamIndex)
		{
		const string &Name = S.m_ParamNames[ParamIndex];
		float Value = ParamValues[ParamIndex];
		if (Name == "gapopen")
			{
			GapOpen = Value;
			printf("gapopen=%8.4g", Value);
			}
		else if (Name == "center")
			{
			Center = Value;
			printf(" center=%8.4g", Value);
			}
		else
			Die("SweeperGetScore bad param '%s'", Name.c_str());
		}

	g_AP.UpdateMx(g_SubstMx_Letter, GapOpen, Center, false);

	g_Bench->Run(g_AP);
	Q = g_Bench->m_MeanQ;
	TC = g_Bench->m_MeanTC;
	printf("  Q=%+6.4f(%6.4f) TC=%+6.4f(%6.4f)",
	  Q - g_TopQ, g_TopQ, TC - g_TopTC, g_TopTC);
	if (S.m_GridCounter != UINT_MAX)
		{
		double Pct = GetPct(S.m_GridCounter, S.m_GridCount);
		printf(" (%.2f%%)", Pct);
		}
	if (S.m_SpatterIter != UINT_MAX)
		printf(" %u[%u/%u]",
		  S.m_SpatterIter, S.m_SpatterTry+1, S.m_SpatterTriesPerIter);
	//double Score = Q + TC;
	double Score = TC;
	if (Score > g_TopScore)
		{
		g_TopScore = Score;
		g_TopQ = Q;
		g_TopTC = TC;
		printf(" TC=%.5f <<", TC);
		printf("\n");
		}
	else
		printf("      \r");
	}

// name,start,delta/name,start,delta...
static void ParseSpatterSpec(const string &Spec,
  vector<string> &Names, vector<float> &Deltas)
	{
	Names.clear();
	Deltas.clear();

	vector<string> Fields;
	Split(Spec, Fields, '/');
	for (uint i = 0; i < SIZE(Fields); ++i)
		{
		const string &Spec1 = Fields[i];
		vector<string> Fields2;
		Split(Spec1, Fields2, ',');
		if (SIZE(Fields2) != 2)
			Die("Bad spec1='%s'", Spec1.c_str());

		const string &Name = Fields2[0];
		float Delta = (float) StrToFloat(Fields2[1]);

		Names.push_back(Name);
		Deltas.push_back(Delta);
		}
	}

void cmd_spatter()
	{
	asserta(optset_warmup_pct);

	Bench B;
	string RefDir = opt(refdir);
	Dirize(RefDir);
	B.Load(opt(spatter), RefDir);

	uint WarmupPct = opt(warmup_pct);
	asserta(WarmupPct > 0 && WarmupPct <= 100);

	if (optset_blosumpct)
		GetSubstMx_Letter_Blosum(opt(blosumpct), g_SubstMx_Letter);
	else if (optset_substmx)
		ReadSubstMx_Letter_FromFile(opt(substmx), g_SubstMx_Letter);
	else
		GetSubstMx_Letter_Blosum(62, g_SubstMx_Letter);

	Bench WarmupB;
	WarmupB.FromSample(B, WarmupPct);

	optset_quiet = true;
	opt_quiet = true;

	g_AP.SetFromCmdLine(false);

	asserta(optset_maxiters);
	asserta(optset_maxfailiters);
	asserta(optset_triesperiter);
	asserta(optset_shrink);
	uint MaxIters = opt(maxiters);				// 16
	uint MaxFailIters = opt(maxfailiters);		// 3
	uint TriesPerIter = opt(triesperiter);		// 32
	float Shrink = (float) opt(shrink);			// 0.6

	asserta(optset_gridspec);
	vector<string> Names;
	vector<float> Goods;
	vector<float> Los;
	vector<float> His;
	vector<uint> Sizes;
	ParseGridSpec(opt(gridspec), Names, Goods, Los, His, Sizes);

	vector<string> Names2;
	vector<float> StartMaxDeltas;
	ParseSpatterSpec(opt(spatterspec), Names2, StartMaxDeltas);
	asserta(SIZE(Names2) == SIZE(Names));
	for (uint i = 0; i < SIZE(Names2); ++i)
		asserta(Names2[i] == Names[i]);

	g_Bench = &WarmupB;
	Sweeper S1;
	S1.m_GetScore = SweeperGetScore;
	S1.SetParamNames(Names);
	S1.SetFev(opt(output1));
	if (!Goods.empty())
		S1.Run1(Goods);
	S1.m_GridNoiseFract = 0.1f;
	S1.ExploreGrid(Los, His, Sizes);

	vector<uint> WarmupIndexes;
	S1.GetDistinctTopIndexes(8, 0.05f, 1.0f, WarmupIndexes);
	const uint WarmupIndexCount = SIZE(WarmupIndexes);
	asserta(WarmupIndexCount > 0);
	vector<vector<float> > StartValueVec;
	for (uint i = 0; i < WarmupIndexCount; ++i)
		{
		uint Index = WarmupIndexes[i];
		const vector<float> &Values = S1.m_ParamValuesVec[Index];
		StartValueVec.push_back(Values);
		}

	ProgressLog("\nWarmup done\n");
	g_TopScore = 0;
	g_Bench = &B;
	Sweeper S2;
	S2.m_GetScore = SweeperGetScore;
	S2.SetParamNames(Names);
	S2.SetFev(opt(fev));
	S2.ExploreSpatter(StartValueVec, StartMaxDeltas,
	  TriesPerIter, MaxIters, MaxFailIters, Shrink);

	S2.LogTop();
	}
