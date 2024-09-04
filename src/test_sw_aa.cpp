#include "muscle.h"
#include "pathscorer.h"
#include "enumpaths.h"
#include "xdpmem.h"

float SWFast_Strings_BLOSUM62(XDPMem &Mem,
  const string &A, const string &B, float Open, float Ext,
  uint &Loi, uint &Loj, uint &Leni, uint &Lenj, string &Path);

static uint g_LA;
static uint g_LB;
static float g_GapOpen = -3;
static float g_GapExt = -1;
static PathScorer_AA_BLOSUM62 *g_PS;
static float g_BestScore;
static string g_BestPath;
static uint g_BestPosA;
static uint g_BestPosB;

static void ClearBrute()
	{
	g_BestScore = -999;
	g_BestPath.clear();
	g_BestPosA = UINT_MAX;
	g_BestPosB = UINT_MAX;
	}

static void OnPath(uint PosA, uint PosB, const string &Path)
	{
	float Score = g_PS->GetLocalScore(PosA, PosB, g_LA, g_LB, Path);
	Log("%10.3g  %5u  %5u  %s\n", Score, PosA, PosB, Path.c_str());
	if (Score > g_BestScore)
		{
		g_BestScore = Score;
		g_BestPath = Path;
		g_BestPosA = PosA;
		g_BestPosB = PosB;
		}
	}

static void Test(const string &A, const string &B)
	{
	ClearBrute();

	g_LA = SIZE(A);
	g_LB = SIZE(B);
	g_PS->m_SeqA = A;
	g_PS->m_SeqB = B;
	EnumPathsLocal(g_LA, g_LB, OnPath);

	XDPMem Mem;
	uint Loi, Loj, Leni, Lenj;
	string SWPath;
	float SWScore = SWFast_Strings_BLOSUM62(Mem, A, B, g_GapOpen,
	  g_GapExt, Loi, Loj, Leni, Lenj, SWPath);
	Log("\n");
	Log("A=%s(%u)", A.c_str(), g_LA);
	Log("  B=%s(%u)", B.c_str(), g_LB);
	Log("  Brute %.3g (%u, %u), SW %.3g (%u, %u) %s %s\n",
	  g_BestScore, g_BestPosA, g_BestPosB,
	  SWScore, Loi, Loj, g_BestPath.c_str(), SWPath.c_str());
	}

static void TestPath(const string &A, const string &B,
  uint PosA, uint PosB, const string &Path)
	{
	ClearBrute();

	g_LA = SIZE(A);
	g_LB = SIZE(B);
	g_PS->m_SeqA = A;
	g_PS->m_SeqB = B;
	g_PS->Trace(true);
	float Score = g_PS->GetLocalScore(PosA, PosB, g_LA, g_LB, Path);
	Log("TestPath %.3g (%u, %u) %s\n",
	  Score, PosA, PosB, Path.c_str());
	}

void cmd_test_sw_aa()
	{
	PathScorer_AA_BLOSUM62 PS;
	PS.m_GapOpen = g_GapOpen;
	PS.m_GapExt = g_GapExt;
	g_PS = &PS;

	Test("SEWWE", "WQW");
	Test("WWESE", "WQW");
	Test("WMW", "WWESE");
	//TestPath("WMW", "WWESE", 0, 0, "MDM");
	//TestPath("WMW", "WWESE", 0, 0, "MIM");
	}
