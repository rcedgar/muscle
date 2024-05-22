#include "muscle.h"
#include "bench.h"
#include "sweeper.h"
#include "m3alnparams.h"

static Bench *g_Bench;
static double g_TopScore = 0;
static double g_TopQ = 0;
static double g_TopTC = 0;
static float g_SubstMx_Letter[20][20];

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

	M3AlnParams AP;
	AP.SetFromCmdLine(false);
	AP.UpdateMx(g_SubstMx_Letter, GapOpen, Center, false);

	if (optset_treeiters)
		AP.m_TreeIters = opt(treeiters);

	g_Bench->Run(AP);
	Q = g_Bench->m_MeanQ;
	TC = g_Bench->m_MeanTC;
	printf("  Q=%+6.4f(%6.4f) TC=%+6.4f(%6.4f)",
	  Q - g_TopQ, g_TopQ, TC - g_TopTC, g_TopTC);
	if (S.m_GridCounter != UINT_MAX)
		{
		double Pct = GetPct(S.m_GridCounter, S.m_GridCount);
		printf(" (%.2f%%)", Pct);
		}
	//double Score = Q + TC;
	double Score = TC;
	if (Score > g_TopScore)
		{
		g_TopScore = Score;
		g_TopQ = Q;
		g_TopTC = TC;
		printf(" <<");
		printf("\n");
		}
	else
		printf("      \r");
	}

// name,good,lo,hi,n/name=good,lo,hi,n...
void ParseGridSpec(const string &Spec,
  vector<string> &Names,
  vector<float> &Goods,
  vector<float> &Los,
  vector<float> &His,
  vector<uint> &Sizes)
	{
	Names.clear();
	Goods.clear();
	Los.clear();
	His.clear();
	Sizes.clear();

	vector<string> Fields;
	Split(Spec, Fields, '/');
	bool DoGoods = true;
	for (uint i = 0; i < SIZE(Fields); ++i)
		{
		const string &Spec1 = Fields[i];
		vector<string> Fields2;
		Split(Spec1, Fields2, ',');
		if (SIZE(Fields2) != 5)
			Die("Bad spec1='%s'", Spec1.c_str());
		const string &Name = Fields2[0];
		if (i == 0 && Fields2[1] == "-")
			DoGoods = false;
		if (DoGoods)
			{
			float Good = (float) StrToFloat(Fields2[1]);
			Goods.push_back(Good);
			}
		float Lo = (float) StrToFloat(Fields2[2]);
		float Hi = (float) StrToFloat(Fields2[3]);
		uint Size = StrToUint(Fields2[4]);
		asserta(Size > 1);
		asserta(Lo != Hi);
		if (Lo < Hi)
			{
			Los.push_back(Lo);
			His.push_back(Hi);
			}
		else
			{
			Los.push_back(Hi);
			His.push_back(Lo);
			}

		Names.push_back(Name);
		Sizes.push_back(Size);
		}
	}

void cmd_sweep()
	{
	asserta(optset_gridspec);
	vector<string> Names;
	vector<float> Goods;
	vector<float> Los;
	vector<float> His;
	vector<uint> Sizes;
	ParseGridSpec(opt(gridspec), Names, Goods, Los, His, Sizes);

	if (optset_blosumpct)
		GetSubstMx_Letter_Blosum(opt(blosumpct), g_SubstMx_Letter);
	else if (optset_substmx)
		ReadSubstMx_Letter_FromFile(opt(substmx), g_SubstMx_Letter);
	else
		GetSubstMx_Letter_Blosum(62, g_SubstMx_Letter);

	Bench B;
	string RefDir = opt(refdir);
	Dirize(RefDir);
	B.Load(g_Arg1, RefDir);
	g_Bench = &B;
	optset_quiet = true;
	opt_quiet = true;

	Sweeper S;
	S.m_GetScore = SweeperGetScore;
	S.SetParamNames(Names);
	S.SetFev(opt(fev));
	if (!Goods.empty())
		S.Run1(Goods);

	S.ExploreGrid(Los, His, Sizes);

	S.LogTop();
	}
