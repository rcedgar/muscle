#if 0
#include "myutils.h"
#include "profpos3.h"
#include "alpha.h"

float ScoreProfPos2(const ProfPos3 &PPA, const ProfPos3 &PPB);

void SetSubstMx2(bool IsNucleo,
  float GapOpen, float GapExt, float Center, float AddCenter,
  uint MulOccs, uint NormAAFreqs, const string &MxName);

static uint GetLetterCounts(const string &Seq,
  vector<uint> &LetterCounts)
	{
	LetterCounts.clear();
	LetterCounts.resize(g_AlphaSize, 0);
	uint GapCount = 0;
	for (uint i = 0; i < SIZE(Seq); ++i)
		{
		char c = Seq[i];
		if (isgap(c))
			{
			++GapCount;
			continue;
			}
		byte Letter = g_CharToLetterAmino[c];
		asserta(Letter < 20);
		LetterCounts[Letter] += 1;
		}
	return GapCount;
	}

static float InnerTest(
  const string &ColA, const string &ColB,
  uint NormAAFreqs, uint MulOccs,
  float Center, float AddCenter,
  uint NA, uint LLA, uint LGA, uint GLA, uint GGA,
  uint NB, uint LLB, uint LGB, uint GLB, uint GGB,
  const vector<uint> &LCA,
  const vector<uint> &LCB
  )
	{
	float GapOpen = -5.5f;
	float GapExt = 0.0f;

	SetSubstMx2(false,
	  GapOpen, GapExt, Center, AddCenter,
	  MulOccs, NormAAFreqs, "B62");

	ProfPos3 PPA;
	ProfPos3 PPB;

	PPA.SetFreqs2(NA, LLA, LGA, GLA, GGA, LCA);
	PPA.SetAAScores();

	PPB.SetFreqs2(NB, LLB, LGB, GLB, GGB, LCB);
	PPB.SetAAScores();

	Log("\n_________________________________\n");
	Log("Norm=%u, MulOccs=%u", NormAAFreqs, MulOccs);
	if (Center != FLT_MAX)
		Log(" Center=%.3g", Center);
	if (AddCenter != FLT_MAX)
		Log(" AddCenter=%.3g", AddCenter);
	Log("\n");

	Log("PPA:\n");
	PPA.LogMe();
	Log("PPB:\n");
	PPB.LogMe();

	float Score = ScoreProfPos2(PPA, PPB);
	Log("Score = %.4g\n", Score);

	Log("&& %7s", ColA.c_str());
	Log(" %7s", ColB.c_str());
	Log(" Norm=%u MulOccs=%u", NormAAFreqs, MulOccs);
	if (Center != FLT_MAX)
		Log("    Center=%6.2f", Center);
	if (AddCenter != FLT_MAX)
		Log(" AddCenter=%6.2f", AddCenter);
	Log(" Score=%8.3f\n", Score);
	
	return Score;
	}

static void Test1(const string &ColA, const string &ColB, float C)
	{
	Log("\n==========================\n");
	Log("ColA = %s\n", ColA.c_str());
	Log("ColB = %s\n", ColB.c_str());

	uint NA = SIZE(ColA);
	uint NB = SIZE(ColB);

	vector<uint> LCA;
	vector<uint> LCB;

	uint LGA = GetLetterCounts(ColA, LCA);
	uint LGB = GetLetterCounts(ColB, LCB);

	uint LLA = NA - LGA;
	uint LLB = NB - LGB;

	uint GLA = 0;
	uint GLB = 0;

	uint GGA = 0;
	uint GGB = 0;

	uint NormAAFreqs = 1;
	uint MulOccs = 1;
	float Center = 0.0f;
	float AddCenter = 0.0f;

	//InnerTest(ColA, ColB, NormAAFreqs, MulOccs,
	//  Center, AddCenter,
	//  NA, LLA, LGA, GLA, GGA,
	//  NB, LLB, LGB, GLB, GGB,
	//  LCA, LCB);

	//NormAAFreqs = 0;
	//MulOccs = 0;
	//Center = 0.0f;
	//AddCenter = 0.0f;

	//InnerTest(ColA, ColB, NormAAFreqs, MulOccs,
	//  Center, AddCenter,
	//  NA, LLA, LGA, GLA, GGA,
	//  NB, LLB, LGB, GLB, GGB,
	//  LCA, LCB);

	NormAAFreqs = 1;
	MulOccs = 1;
	Center = C;
	AddCenter = 0;
	float Score1 = InnerTest(ColA, ColB, NormAAFreqs, MulOccs,
	  Center, AddCenter,
	  NA, LLA, LGA, GLA, GGA,
	  NB, LLB, LGB, GLB, GGB,
	  LCA, LCB);

	NormAAFreqs = 0;
	MulOccs = 0;
	Center = 0;
	AddCenter = C;
	float Score2 = InnerTest(ColA, ColB, NormAAFreqs, MulOccs,
	  Center, AddCenter,
	  NA, LLA, LGA, GLA, GGA,
	  NB, LLB, LGB, GLB, GGB,
	  LCA, LCB);

	if (feq(Score1, Score2))
		Log(" == Equivalent OK ==\n");
	else
		Die("Score1 %.5g, Score2 %.5g\n", Score1, Score2);
	}

void cmd_scoretest()
	{
	opt(scoretest);

	SetAlpha(ALPHA_Amino);
	//extern float g_CENTER;
	//g_CENTER = 0;
	Test1("C-", "CC", 0.5);
	Test1("SEQ", "SEQ", 1.0);
	Test1("AAC-", "MCC", 2.0);
	Test1("AAC-", "SEQ---", 2.5);
	}
#endif // 0
