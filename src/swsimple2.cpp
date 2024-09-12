#include "myutils.h"
#include "mx.h"
#include "pathscorer.h"
#include "allocmx.h"
#include "xdpmem.h"
#include "tracebit.h"

inline bool myfeq(double x, double y)
    {
	if (x == y)
		return true;
	if (x < -999 && y < -999)
		return true;
    double X = fabs(x);
    double Y = fabs(y);
    double Max = max(X, Y);
    double Diff = fabs(X-Y);
    return Diff < Max*0.01;
    }

static uint trace_i = UINT_MAX;
static uint trace_j = UINT_MAX;
static char trace_mx = 0;

float SWSimpleFwdMDI(PathScorer &PS, uint &LoA, uint &LoB, string &Path,
  vector<vector<float> > &FwdM,
  vector<vector<float> > &FwdD,
  vector<vector<float> > &FwdI,
  vector<vector<char> > &TBM,
  vector<vector<char> > &TBD,
  vector<vector<char> > &TBI);

static void Logx(float x)
	{
	if (x == MINUS_INFINITY)
		Log("  %7.7s", "*");
	else if (x == FLT_MAX)
		Log("  %7.7s", "&");
	else if (x == UNINIT)
		Log("  %7.7s", " ");
	else
		Log("  %7.2f", x);
	}

static void LogMx(const char *Name, const vector<vector<float> > &Mx)
	{
	const uint LA = SIZE(Mx);
	const uint LB = SIZE(Mx[0]);
	Log("\nLogMx(%s)\n", Name);
	Log("      ");
	for (uint j = 0; j < LB; ++j)
		Log("  %7u", j);
	Log("\n");
	for (uint i = 0; i < LA; ++i)
		{
		Log("%3u | ", i);
		for (uint j = 0; j < LB; ++j)
			Logx(Mx[i][j]);
		Log("\n");
		}
	Log("\n");
	}

static void LogTBMx(const char *Name, const vector<vector<char> > &Mx)
	{
	const uint LA = SIZE(Mx);
	const uint LB = SIZE(Mx[0]);
	Log("\nLogMx(%s)\n", Name);
	Log("      ");
	for (uint j = 0; j < LB; ++j)
		Log("  %u", j%10);
	Log("\n");
	for (uint i = 0; i < LA; ++i)
		{
		Log("%3u | ", i);
		for (uint j = 0; j < LB; ++j)
			Log("%c", Mx[i][j]);
		Log("\n");
		}
	Log("\n");
	}

static void CmpMx(char c, const vector<vector<float> > &SM,
  const vector<vector<float> > &M)
 	{
	uint LA = SIZE(SM);
	asserta(SIZE(M) == LA);
	uint LB = SIZE(SM[0]);
	asserta(SIZE(M[0]) == LB);

	if (c != 'M')
		{
		--LA;
		--LB;
		}
	for (uint i = 1; i < LA; ++i)
		{
		for (uint j = 1; j < LB; ++j)
			{
			float sx = SM[i][j];
			float x = M[i][j];
			if (!myfeq(sx, x))
				{
				Log("%c i=%u j=%u %.3f %.3f\n", c, i, j, sx, x);
				Die("CmpMx");
				}
			}
		}
	}

#define CMPM(i, j)	\
	{ \
    float x = FwdM[i+1][j+1]; \
	float sx = Simple_FwdM[i+1][j+1]; \
	if (!myfeq(x, sx)) \
		{ \
		LogMx("Simple_M", Simple_FwdM); \
		LogMx("Fast_M", FwdM); \
		LogMx("Simple_D", Simple_FwdD); \
		LogMx("Fast_D", FwdD); \
		LogMx("Simple_I", Simple_FwdI); \
		LogMx("Fast_I", FwdI); \
		Die("CMP i+1=%u j+1=%u FwdM=%.3f SimpleM=%.3f", i+1, j+1, x, sx); \
		} \
	}

#define CMPD(i, j)	\
	{ \
    float x = FwdD[i+1][j+1]; \
	float sx = Simple_FwdD[i+1][j+1]; \
	if (!myfeq(x, sx) && i > 0 && j > 0 && i+1 < LA && j+1 < LB) \
		{ \
		LogMx("Simple_M", Simple_FwdM); \
		LogMx("Fast_M", FwdM); \
		LogMx("Simple_D", Simple_FwdD); \
		LogMx("Fast_D", FwdD); \
		LogMx("Simple_I", Simple_FwdI); \
		LogMx("Fast_I", FwdI); \
		Die("CMP i+1=%u j+1=%u FwdD=%.3f SimpleD=%.3f", i+1, j+1, x, sx); \
		} \
	}

#define CMPI(i, j)	\
	{ \
    float x = FwdI[i+1][j+1]; \
	float sx = Simple_FwdI[i+1][j+1]; \
	if (!myfeq(x, sx) && i > 0 && j > 0 && i+1 < LA && j+1 < LB) \
		{ \
		LogMx("Simple_M", Simple_FwdM); \
		LogMx("Fast_M", FwdM); \
		LogMx("Simple_D", Simple_FwdD); \
		LogMx("Fast_D", FwdD); \
		LogMx("Simple_I", Simple_FwdI); \
		LogMx("Fast_I", FwdI); \
		Die("CMP i+1=%u j+1=%u FwdI=%.3f SimpleI=%.3f", i+1, j+1, x, sx); \
		} \
	}

float SWSimple2(XDPMem &Mem, PathScorer &PS, uint &LoA, uint &LoB, string &Path)
	{
	vector<vector<float> > Simple_FwdM;
	vector<vector<float> > Simple_FwdD;
	vector<vector<float> > Simple_FwdI;
	vector<vector<char> > Simple_TBM;
	vector<vector<char> > Simple_TBD;
	vector<vector<char> > Simple_TBI;
	string Simple_Path;
	uint Simple_LoA;
	uint Simple_LoB;
	float Simple_Score =
	  SWSimpleFwdMDI(PS, Simple_LoA, Simple_LoB, Simple_Path,
		Simple_FwdM, Simple_FwdD, Simple_FwdI,
		Simple_TBM, Simple_TBD, Simple_TBI);
	LogMx("Simple_M", Simple_FwdM);
	LogMx("Simple_D", Simple_FwdD);
	LogMx("Simple_I", Simple_FwdI);
	LogTBMx("Simple_TBM", Simple_TBM);
	LogTBMx("Simple_TBD", Simple_TBD);
	LogTBMx("Simple_TBI", Simple_TBI);

	const uint LA = PS.GetLA();
	const uint LB = PS.GetLB();

	float SWFast_SMx(XDPMem &Mem, const Mx<float> &SMx,
	  float Open, float Ext, uint &Loi, uint &Loj, uint &Leni, uint &Lenj,
	  string &Path);
	Mx<float> SMx;
	SMx.Alloc(LA, LB);
	for (uint i = 0; i < LA; ++i)
		{
		for (uint j = 0; j < LB; ++j)
			{
			float x = PS.GetMatchScore(i, j);
			SMx.Put(i, j, x);
			}
		}
	float Open = PS.GetScoreMD(0, 0);
	float Ext = PS.GetScoreDD(0, 0);
	uint SMx_Loi, SMx_Loj, SMx_Leni, SMx_Lenj;
	string SMx_Path;
	float SMx_Score = SWFast_SMx(Mem, SMx, Open, Ext,
	  SMx_Loi, SMx_Loj, SMx_Leni, SMx_Lenj, SMx_Path);

	vector<vector<float> > FwdM;
	vector<vector<float> > FwdD;
	vector<vector<float> > FwdI;
	AllocMx<float>(FwdM, LA+1, LB+1, UNINIT);
	AllocMx<float>(FwdD, LA+1, LB+1, UNINIT);
	AllocMx<float>(FwdI, LA+1, LB+1, UNINIT);

	Mem.Alloc(LA+32, LB+32);

	float *Mrow = Mem.GetDPRow1();
	float *Drow = Mem.GetDPRow2();
	byte **TB = Mem.GetTBBit();

// Use Mrow[-1], so...
	Mrow[-1] = MINUS_INFINITY;

	for (uint j = 0; j <= LB; ++j)
		{
		Mrow[j] = MINUS_INFINITY;
		Drow[j] = MINUS_INFINITY;
		}
	
	float BestScore = 0;
	uint Besti = UINT_MAX;
	uint Bestj = UINT_MAX;

#define eq(x, y)	{ if (!myfeq((x), (y))) \
  { LogMx("FwdM", FwdM); LogMx("FwdD", FwdD); LogMx("FwdI", FwdI); \
	Die("%d: i=%u j=%u myfeq(%.3g, %.3g) %s %s", __LINE__, i, j, (x), (y), #x, #y);} }

	for (uint i = 0; i < LA; ++i)
		{
		float Is = MINUS_INFINITY;
		byte *TBrow = TB[i];
		float PrevSavedM = MINUS_INFINITY;
		for (uint j = 0; j < LB; ++j)
			{
			byte TraceBits = 0;
			float SavedM = Mrow[j+1];
			eq(PrevSavedM, Simple_FwdM[i][j]);
			eq(SavedM, Simple_FwdM[i][j+1]);

			////////////////////////////
			//	M[i+1][j+1] = max
			//		M[i][j] + tMM + sub,
			//		D[i][j] + tDM + sub,
			//		I[i][j] + tIM + sub
			////////////////////////////

			// Match input:
			//		PrevSavedM = M[i][j]
			//		Drow[j] = D[i][j]
			//		Is = I[i][j]
			eq(PrevSavedM, Simple_FwdM[i][j]);
			eq(Drow[j], Simple_FwdD[i][j]);
			eq(Is, Simple_FwdI[i][j]);

			float mm = PrevSavedM + PS.GetScoreMM(i, j);
			float dm = Drow[j] + PS.GetScoreDM(i, j);
			float im = Is + PS.GetScoreIM(i, j);

			float t = mm;
			if (dm > t)
				{
				t = dm;
				TraceBits = TRACEBITS_DM;
				}
			if (im > t)
				{
				t = im;
				TraceBits = TRACEBITS_IM;
				}
			if (t < 0)
				{
				TraceBits = TRACEBITS_SM;
				t = 0;
				}
			float m = t + PS.GetMatchScore(i, j);
			if (m > BestScore)
				{
				BestScore = m;
				Besti = i;
				Bestj = j;
				}

			// Match output (I = i+1, i will have been incremented when used):
			//		Mrow[j+1] = M[I][j+1]
			Mrow[j+1] = m;
			eq(Mrow[j+1], Simple_FwdM[i+1][j+1]);
			FwdM[i+1][j+1] = m;

			////////////////////////////
			//	D[i+1][j] = max
			//		M[i][j] + tMD,
			//		D[i][j] + tDD
			////////////////////////////

			// Delete input:
			//	PrevSavedM = M[i][j]
			//	Drow[j] = D[i][j]
			eq(PrevSavedM, Simple_FwdM[i][j]);
			eq(Drow[j], Simple_FwdD[i][j]);
			float md = PrevSavedM + PS.GetScoreMD(i, j);
			float dd = Drow[j] + PS.GetScoreDD(i, j);

			// Delete output:
			//		Drow[j] = D[I][j]
			//float d = max(md, dd);
			float d = dd;
			if (md > dd)
				{
				d = md;
				TraceBits |= TRACEBITS_MD;
				}
			Drow[j] = d;
			eq(Drow[j], Simple_FwdD[i+1][j]);
			FwdD[i+1][j] = d;

			// Insert input:
			//	PrevSavedM = M[i][j]
			//	Is = I[i][j]
			eq(PrevSavedM, Simple_FwdM[i][j]);
			eq(Is, Simple_FwdI[i][j]);
			float mi = PrevSavedM + PS.GetScoreMI(i, j);
			float ii = Is + PS.GetScoreII(i, j);
			//Is = max(mi, ii);
			Is = ii;
			if (mi > ii)
				{
				Is = mi;
				TraceBits |= TRACEBITS_MI;
				}
			// Insert output:
			//	Is = I[i][j+1]
			eq(Is, Simple_FwdI[i][j+1]);
			FwdI[i][j+1] = Is;

			PrevSavedM = SavedM;
			TBrow[j] = TraceBits;
			}
		//Log("i=%u Mrow=\n", i);
		//for (uint j = 0; j <= LB; ++j)
		//	Logx(Mrow[j]);
		//Log("\n");
		//Log("i=%u Drow=\n", i);
		//for (uint j = 0; j <= LB; ++j)
		//	Logx(Drow[j]);
		//Log("\n");
		}

	asserta(myfeq(BestScore, Simple_Score));
	asserta(myfeq(BestScore, SMx_Score));

//////////////////////////////////////////////////////////////
	if (BestScore <= 0.0f)
		return 0.0f;

	LogTBSW("Simple2", Mem, LA, LB);
	uint Leni, Lenj;
	TraceBackBitSW(Mem, LA, LB, Besti, Bestj,
	  Leni, Lenj, Path);
	asserta(Besti+1 >= Leni);
	asserta(Bestj+1 >= Lenj);
	asserta(Path == Simple_Path);

	LoA = Besti + 1 - Leni;
	LoB = Bestj + 1 - Lenj;

	//LogMx("Simple_M", Simple_FwdM);
	//LogMx("Fast_M", FwdM);

	//LogMx("Simple_D", Simple_FwdD);
	//LogMx("Fast_D", FwdD);

	//LogMx("Simple_I", Simple_FwdI);
	//LogMx("Fast_I", FwdI);

	//CmpMx('M', Simple_FwdM, FwdM);
	//CmpMx('D', Simple_FwdD, FwdD);
	//CmpMx('I', Simple_FwdI, FwdI);

	return BestScore;
	}

void cmd_swsimple2()
	{
	string A = "SEQVENCE";
	string B = "QVEN";

	PathScorer_AA_BLOSUM62 PS;
	PS.m_GapOpen = -2.2f;
	PS.m_GapExt = -0.5f;
	PS.m_SeqA = A;
	PS.m_SeqB = B;
	PS.m_LA = SIZE(A);
	PS.m_LB = SIZE(B);

	XDPMem Mem;
	uint LoA, LoB;
	string Path;
	float Score = SWSimple2(Mem, PS, LoA, LoB, Path);
	Log("Score %.3g path=%s\n", Score, Path.c_str());
	}
