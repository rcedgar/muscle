#include "myutils.h"
#include "xdpmem.h"
#include "tracebit.h"
#include "objmgr.h"
#include "pathinfo.h"
#include "mega.h"

/***
https://drive5workspace.slack.com/archives/D076HKBFAM6/p1725457191133919
gapopen=0.847836 gapext=0.105778 termgapopen=0.000000 termgapext=0.097557
TC = [0.8877]

-gapopen=0.85 -gapext=0.10 -termgapopen=0.0 -termgapext=0.10
TC = [0.91]
***/

void TraceBackBitMem(XDPMem &Mem, unsigned LA, unsigned LB, char State, PathInfo &PI);

static float g_AP_LOpenA;
static float g_AP_LOpenB;
static float g_AP_LExtA;
static float g_AP_LExtB;

static float g_AP_ROpenA;
static float g_AP_ROpenB;
static float g_AP_RExtA;
static float g_AP_RExtB;

static float g_AP_OpenA;
static float g_AP_OpenB;
static float g_AP_ExtA;
static float g_AP_ExtB;

static void SetGaps(float IntOpen, float IntExt, float TermOpen, float TermExt)
	{
	g_AP_LOpenA = TermOpen;
	g_AP_LOpenB = TermOpen;
	g_AP_LExtA = TermExt;
	g_AP_LExtB = TermExt;

	g_AP_ROpenA = TermOpen;
	g_AP_ROpenB = TermOpen;
	g_AP_RExtA = TermExt;
	g_AP_RExtB = TermExt;

	g_AP_OpenA = IntOpen;
	g_AP_OpenB = IntOpen;
	g_AP_ExtA = IntExt;
	g_AP_ExtB = IntExt;
	}

static float ViterbiMega(XDPMem &Mem, const vector<vector<byte> > &ProfA,
   const vector<vector<byte> > &ProfB, PathInfo &PI)
	{
	const uint LA = SIZE(ProfA);
	const uint LB = SIZE(ProfB);
	if (LA*LB > 100*1000*1000)
		Die("ViterbiMega, seqs too long LA=%u, LB=%u", LA, LB);

	Mem.Alloc(LA, LB);
	PI.Alloc2(LA, LB);
	

	float OpenA = g_AP_LOpenA;
	float ExtA = g_AP_LExtA;

	float *Mrow = Mem.GetDPRow1();
	float *Drow = Mem.GetDPRow2();
	byte **TB = Mem.GetTBBit();

// Use Mrow[-1], so...
	Mrow[-1] = MINUS_INFINITY;
	for (unsigned j = 0; j <= LB; ++j)
		{
		Mrow[j] = MINUS_INFINITY;
		Drow[j] = MINUS_INFINITY;
		}
	
// Main loop
	float M0 = float (0);
	for (unsigned i = 0; i < LA; ++i)
		{
		//byte a = A[i];
		//const float *MxRow = Mx[a];
		float OpenB = g_AP_LOpenB;
		float ExtB = g_AP_LExtB;
		float I0 = MINUS_INFINITY;

		byte *TBrow = TB[i];
		for (unsigned j = 0; j < LB; ++j)
			{
			//byte b = B[j];
			byte TraceBits = 0;
			float SavedM0 = M0;

		// MATCH
			{
		// M0 = DPM[i][j]
		// I0 = DPI[i][j]
		// Drow[j] = DPD[i][j]

			float xM = M0;
			if (Drow[j] > xM)
				{
				xM = Drow[j];
				TraceBits = TRACEBITS_DM;
				}
			if (I0 > xM)
				{
				xM = I0;
				TraceBits = TRACEBITS_IM;
				}
			M0 = Mrow[j];

			float MatchScore = Mega::GetMatchScore_LogOdds(ProfA, i, ProfB, j);
			Mrow[j] = xM + MatchScore;
		// Mrow[j] = DPM[i+1][j+1])
			}
			
		// DELETE
			{
		// SavedM0 = DPM[i][j]
		// Drow[j] = DPD[i][j]
			float md = SavedM0 + OpenB;
			Drow[j] += ExtB;
			if (md >= Drow[j])
				{
				Drow[j] = md;
				TraceBits |= TRACEBITS_MD;
				}
		// Drow[j] = DPD[i+1][j]
			}
			
		// INSERT
			{
		// SavedM0 = DPM[i][j]
		// I0 = DPI[i][j]
			float mi = SavedM0 + OpenA;
			I0 += ExtA;
			if (mi >= I0)
				{
				I0 = mi;
				TraceBits |= TRACEBITS_MI;
				}
		// I0 = DPI[i][j+1]
			}
			
			OpenB = g_AP_OpenB;
			ExtB = g_AP_ExtB;
			
			TBrow[j] = TraceBits;
			}
		
	// Special case for end of Drow[]
		{
	// M0 = DPM[i][LB]
	// Drow[LB] = DPD[i][LB]
		
		TBrow[LB] = 0;
		float md = M0 + g_AP_ROpenB;
		Drow[LB] += g_AP_RExtB;
		if (md >= Drow[LB])
			{
			Drow[LB] = md;
			TBrow[LB] = TRACEBITS_MD;
			}
	// Drow[LB] = DPD[i+1][LB]
		}
		
		M0 = MINUS_INFINITY;

		OpenA = g_AP_OpenA;
		ExtA = g_AP_ExtA;
		}
	
// Special case for last row of DPI
	byte *TBrow = TB[LA];
	float I1 = MINUS_INFINITY;
	for (unsigned j = 1; j < LB; ++j)
		{
	// Mrow[j-1] = DPM[LA][j]
	// I1 = DPI[LA][j]
		
		TBrow[j] = 0;
		float mi = Mrow[int(j)-1] + g_AP_ROpenA;
		I1 += g_AP_RExtA;
		if (mi > I1)
			{
			I1 = mi;
			TBrow[j] = TRACEBITS_MI;
			}
		}
	
	float FinalM = Mrow[LB-1];
	float FinalD = Drow[LB];
	float FinalI = I1;
// FinalM = DPM[LA][LB]
// FinalD = DPD[LA][LB]
// FinalI = DPI[LA][LB]
	
	float Score = FinalM;
	byte State = 'M';
	if (FinalD > Score)
		{
		Score = FinalD;
		State = 'D';
		}
	if (FinalI > Score)
		{
		Score = FinalI;
		State = 'I';
		}

	TraceBackBitMem(Mem, LA, LB, State, PI);

	return Score;
	}

void AlignMega2()
	{
// pymoo optimization
// gapopen=0.847836 gapext=0.105778 termgapopen=0.000000 termgapext=0.097557
	float IntOpen = float(-optd(gapopen, 0.85));
	float IntExt = float(-optd(gapext, 0.10));
	float TermOpen = float(-optd(termgapopen, 0.0));
	float TermExt = float(-optd(termgapext, 0.10));
	SetGaps(IntOpen, IntExt, TermOpen, TermExt);

	const uint ProfileCount = Mega::GetProfileCount();
	if (ProfileCount != 2)
		Die("AlignMega2(): %u structures found, 2 required", ProfileCount);

	XDPMem Mem;
	PathInfo *PI = ObjMgr::GetPathInfo();
	const vector<vector<byte> > &ProfA = Mega::GetProfile(0);
	const vector<vector<byte> > &ProfB = Mega::GetProfile(1);
	float Score = ViterbiMega(Mem, ProfA, ProfB, *PI);
	if (!optset_output)
		return;

	FILE *fOut = CreateStdioFile(opt(output));

	const string &SeqA = Mega::m_Seqs[0];
	const string &SeqB = Mega::m_Seqs[1];

	const uint ColCount = PI->GetColCount();
	const char *Path = PI->GetPath();
	uint PosA = 0;
	uint PosB = 0;

	fprintf(fOut, ">%s\n", Mega::m_Labels[0].c_str());
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];
		if (c == 'M' || c == 'D')
			fprintf(fOut, "%c", SeqA[PosA++]);
		else
			fprintf(fOut, "-");
		}
	fprintf(fOut, "\n");
	fprintf(fOut, ">%s\n", Mega::m_Labels[1].c_str());
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];
		if (c == 'M' || c == 'I')
			fprintf(fOut, "%c", SeqB[PosB++]);
		else
			fprintf(fOut, "-");
		}
	fprintf(fOut, "\n");

	CloseStdioFile(fOut);
	}

void cmd_mega2()
	{
	Mega::FromFile(g_Arg1);
	AlignMega2();
	}
