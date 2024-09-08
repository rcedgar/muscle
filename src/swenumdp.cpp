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
static float g_BestScore;

static void OnPath(uint PosA, uint PosB, const string &Path)
	{
	float Score = g_PS->GetLocalScore(PosA, PosB, Path);
	uint NA = GetNA(Path);
	uint NB = GetNA(Path);
	uint IxA = PosA + NA;
	uint IxB = PosA + NB;
	vector<vector<float> > &FwdM = *g_FwdM;
	asserta(IxA < SIZE(FwdM));
	vector<float> &RowA = FwdM[IxA];
	asserta(IxB < SIZE(RowA));
	if (RowA[IxB] == FLT_MAX || Score > RowA[IxB])
		RowA[IxB] = Score;
	}

float SWEnumDP(PathScorer &PS, uint LA, uint LB,
  vector<vector<float> > &FwdM)
	{
	g_PS = &PS;
	g_FwdM = &FwdM;
	g_LA = LA;
	g_LB = LB;
	g_BestScore = 0;
	AllocMx(FwdM, LA+1, LB+1, FLT_MAX);
	EnumPathsLocal(LA, LB, OnPath);
	return g_BestScore;
	}
