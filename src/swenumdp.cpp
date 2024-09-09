#include "myutils.h"
#include "enumpaths.h"
#include "pathscorer.h"
#include "allocmx.h"

uint GetNA(const string &Path);
uint GetNB(const string &Path);

static vector<vector<float> > *g_FwdM;
static PathScorer *g_PS;
static uint g_LA;
static uint g_LB;
static uint g_BestPosA;
static uint g_BestPosB;
static float g_BestScore;
static string g_BestPath;

static void OnPath(uint PosA, uint PosB, const string &Path)
	{
	float Score = g_PS->GetLocalScore(PosA, PosB, Path);
	if (Score > g_BestScore)
		{
		g_BestScore = Score;
		g_BestPosA = PosA;
		g_BestPosB = PosB;
		g_BestPath = Path;
		}
	uint NA = GetNA(Path);
	uint NB = GetNB(Path);
	uint IxA = PosA + NA;
	uint IxB = PosB + NB;
	vector<vector<float> > &FwdM = *g_FwdM;
	asserta(IxA < SIZE(FwdM));
	vector<float> &RowA = FwdM[IxA];
	asserta(IxB < SIZE(RowA));
	if (RowA[IxB] == FLT_MAX || Score > RowA[IxB])
		RowA[IxB] = Score;
	}

float SWEnumDPFwdM(PathScorer &PS, uint LA, uint LB, uint &LoA, uint &LoB,
  string &Path, vector<vector<float> > &FwdM)
	{
	g_BestScore = 0;
	g_BestPosA = UINT_MAX;
	g_BestPosB = UINT_MAX;
	g_PS = &PS;
	g_FwdM = &FwdM;
	g_LA = LA;
	g_LB = LB;
	g_BestScore = 0;
	AllocMx(FwdM, LA+1, LB+1, FLT_MAX);
	for (uint i = 0; i <= LA; ++i)
		FwdM[i][0] = 0;
	for (uint j = 0; j <= LB; ++j)
		FwdM[0][j] = 0;
	EnumPathsLocal(LA, LB, OnPath);
	LoA = g_BestPosA;
	LoB = g_BestPosB;
	Path = g_BestPath;
	return g_BestScore;
	}

float SWEnumDP(PathScorer &PS, uint LA, uint LB, uint &LoA, uint &LoB,
  string &Path)
	{
	vector<vector<float> > FwdM;
	float Score = SWEnumDPFwdM(PS, LA, LB, LoA, LoB, Path, FwdM);
	return Score;
	}
