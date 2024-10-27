#include "myutils.h"
#include "mx.h"
#include "pathscorer.h"
#include "allocmx.h"

float SWSimpleFwdMDI(PathScorer &PS, uint &LoA, uint &LoB, string &Path,
  vector<vector<float> > &FwdM,
  vector<vector<float> > &FwdD,
  vector<vector<float> > &FwdI,
  vector<vector<char> > &TBM,
  vector<vector<char> > &TBD,
  vector<vector<char> > &TBI)
	{
	uint LA = PS.GetLA();
	uint LB = PS.GetLB();

	AllocMx(FwdM, LA+1, LB+1, FLT_MAX);
	AllocMx(FwdD, LA+1, LB+1, FLT_MAX);
	AllocMx(FwdI, LA+1, LB+1, FLT_MAX);

	AllocMx(TBM, LA+1, LB+1, '?');
	AllocMx(TBD, LA+1, LB+1, '?');
	AllocMx(TBI, LA+1, LB+1, '?');

	for (uint i = 0; i <= LA; ++i)
		{
		FwdM[i][0] = MINUS_INFINITY;
		FwdD[i][0] = MINUS_INFINITY;
		FwdI[i][0] = MINUS_INFINITY;
		TBM[i][0] = 'S';
		TBD[i][0] = '?';
		TBI[i][0] = '?';
		}

	for (uint j = 0; j <= LB; ++j)
		{
		FwdM[0][j] = MINUS_INFINITY;
		FwdD[0][j] = MINUS_INFINITY;
		FwdI[0][j] = MINUS_INFINITY;
		TBM[0][j] = 'S';
		TBD[0][j] = '?';
		TBI[0][j] = '?';
		}

// Main loop
	float BestScore = 0;
	uint Besti = UINT_MAX;
	uint Bestj = UINT_MAX;
	for (uint i = 0; i < LA; ++i)
		{
		//byte a = A[i];
		//const float *SubstMxRow = g_SubstMx[a];
		for (uint j = 0; j < LB; ++j)
			{
		// xM
			{
			float m = PS.GetMatchScore(i, j);
			float SM = m;
			float MM = FwdM[i][j] + PS.GetScoreMM(i, j) + m;
			float DM = FwdD[i][j] + PS.GetScoreDM(i, j) + m;
			float IM = FwdI[i][j] + PS.GetScoreIM(i, j) + m;

			float s = MM;
			char t = 'M';
			if (DM > s)
				{
				s = DM;
				t = 'D';
				}
			if (IM > s)
				{
				s = IM;
				t = 'I';
				}
			if (SM >= s)
				{
				s = SM;
				t = 'S';
				}

			FwdM[i+1][j+1] = s;
			TBM[i+1][j+1] = t;
			if (s > BestScore)
				{
				BestScore = s;
				Besti = i+1;
				Bestj = j+1;
				}
			}
			
		// xD
			{
			float MD = FwdM[i][j+1] + PS.GetScoreMD(i, j+1);
			float DD = FwdD[i][j+1] + PS.GetScoreDD(i, j+1);
			if (MD >= DD)
				{
				FwdD[i+1][j+1] = MD;
				TBD[i+1][j+1] = 'M';
				}
			else
				{
				FwdD[i+1][j+1] = DD;
				TBD[i+1][j+1] = 'D';
				}
			}
			
		// xI
			{
			float MI = FwdM[i+1][j] + PS.GetScoreMI(i+1, j);
			float II = FwdI[i+1][j] + PS.GetScoreII(i+1, j);
			if (MI >= II)
				{
				FwdI[i+1][j+1] = MI;
				TBI[i+1][j+1] = 'M';
				}
			else
				{
				FwdI[i+1][j+1] = II;
				TBI[i+1][j+1] = 'I';
				}
			}
			}
		}

	if (Besti == UINT_MAX)
		return 0;

	uint i = Besti;
	uint j = Bestj;

	Path.clear();
	char State = 'M';
	for (;;)
		{
		Path += State;
		switch (State)
			{
		case 'M':
			{
			State = TBM[i][j];
			--i;
			--j;
			break;
			}
		case 'D':
			{
			State = TBD[i][j];
			--i;
			break;
			}
		case 'I':
			{
			State = TBI[i][j];
			--j;
			break;
			}
		default:
			asserta(false);
			}

		if (i == 0 || j == 0 || State == 'S')
			break;
		}

	reverse(Path.begin(), Path.end());

	LoA = i;
	LoB = j;

	asserta(Besti >= LoA);
	asserta(Bestj >= LoB);

	return BestScore;
	}

float SWSimpleFwdM(PathScorer &PS, uint &LoA, uint &LoB, string &Path,
  vector<vector<float> > &FwdM)
	{
	vector<vector<float> > FwdD;
	vector<vector<float> > FwdI;
	vector<vector<char> > TBM;
	vector<vector<char> > TBD;
	vector<vector<char> > TBI;
	float Score = SWSimpleFwdMDI(PS, LoA, LoB, Path, FwdM, FwdD, FwdI,
	  TBM, TBD, TBI);
	return Score;
	}

float SWSimple(PathScorer &PS, uint &LoA, uint &LoB, string &Path)
	{
	vector<vector<float> > FwdM;
	vector<vector<float> > FwdD;
	vector<vector<float> > FwdI;
	vector<vector<char> > TBM;
	vector<vector<char> > TBD;
	vector<vector<char> > TBI;
	float Score = SWSimpleFwdMDI(PS, LoA, LoB, Path, FwdM, FwdD, FwdI,
	  TBM, TBD, TBI);
	return Score;
	}
