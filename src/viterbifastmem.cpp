#include "myutils.h"
#include "xdpmem.h"
#include "tracebit.h"
#include "pathinfo.h"
#include "m3alnparams.h"
#include "viterbiparams.h"
#include "alpha.h"
#include <mutex>

float g_Viterbi_GapOpen = FLT_MAX;
float g_Viterbi_GapExt = FLT_MAX;
float g_Viterbi_TermGapOpen = FLT_MAX;
float g_Viterbi_TermGapExt = FLT_MAX;
float **g_Viterbi_SubstMx_Char = 0;
static float Blosum62[20][20];

void GetSubstMx_Letter_Blosum(uint PctId, float Mx[20][20]);

static bool g_InitParamsDone = false;
static void SetParams()
	{
	static mutex mut;
	mut.lock();
	if (!g_InitParamsDone)
		{
		g_Viterbi_GapOpen = -6;
		g_Viterbi_GapExt = -1;
		g_Viterbi_TermGapExt = -1;
		g_Viterbi_TermGapOpen = -6;
		g_InitParamsDone = true;
		GetSubstMx_Letter_Blosum(62, Blosum62);
		g_Viterbi_SubstMx_Char = myalloc(float *, 256);
		for (uint i = 0; i < 256; ++i)
			{
			g_Viterbi_SubstMx_Char[i] = myalloc(float, 256);
			for (uint j = 0; j < 256; ++j)
				g_Viterbi_SubstMx_Char[i][j] = 0;
			}

		for (uint i = 0; i < 20; ++i)
			{
			char ci = g_LetterToCharAmino[i];
			for (uint j = 0; j < 20; ++j)
				{
				char cj = g_LetterToCharAmino[i];
				float Score = Blosum62[i][j];
				g_Viterbi_SubstMx_Char[ci][cj] = Score;
				}
			}
		g_InitParamsDone = true;
		}
	mut.unlock();
	}

float ViterbiFastMem(XDPMem &Mem, const byte *A, uint LA,
  const byte *B, uint LB, PathInfo &PI)
	{
	if (LA*LB > 100*1000*1000)
		Die("ViterbiFastMem, seqs too long LA=%u, LB=%u", LA, LB);
	SetParams();

	Mem.Alloc(LA, LB);
	PI.Alloc2(LA, LB);
	
	float OpenA = g_Viterbi_TermGapOpen;
	float ExtA = g_Viterbi_TermGapExt;

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
		byte a = A[i];
		const float *MxRow = g_Viterbi_SubstMx_Char[a];
		float OpenB = g_Viterbi_TermGapOpen;
		float ExtB = g_Viterbi_TermGapExt;
		float I0 = MINUS_INFINITY;

		byte *TBrow = TB[i];
		for (unsigned j = 0; j < LB; ++j)
			{
			byte b = B[j];
			byte TraceBits = 0;
			float SavedM0 = M0;

		// MATCH
			{
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

		// Mrow[j] = xM + (toupper(a) == toupper(b) ? MATCH_SCORE : MISMATCH_SCORE);
			Mrow[j] = xM + MxRow[b];
			}
			
		// DELETE
			{
			float md = SavedM0 + OpenB;
			Drow[j] += ExtB;
			if (md >= Drow[j])
				{
				Drow[j] = md;
				TraceBits |= TRACEBITS_MD;
				}
			}
			
		// INSERT
			{
			float mi = SavedM0 + OpenA;
			I0 += ExtA;
			if (mi >= I0)
				{
				I0 = mi;
				TraceBits |= TRACEBITS_MI;
				}
		// I0 = DPI[i][j+1]
			}
			
			OpenB = g_Viterbi_GapOpen;
			ExtB = g_Viterbi_GapExt;
			
			TBrow[j] = TraceBits;
			}
		
	// Special case for end of Drow[]
		{
		TBrow[LB] = 0;
		float md = M0 + g_Viterbi_TermGapOpen;
		Drow[LB] += g_Viterbi_TermGapExt;
		if (md >= Drow[LB])
			{
			Drow[LB] = md;
			TBrow[LB] = TRACEBITS_MD;
			}
		}
		
		M0 = MINUS_INFINITY;

		OpenA = g_Viterbi_GapOpen;
		ExtA = g_Viterbi_GapExt;
		}
	
// Special case for last row of DPI
	byte *TBrow = TB[LA];
	float I1 = MINUS_INFINITY;
	for (unsigned j = 1; j < LB; ++j)
		{
	// Mrow[j-1] = DPM[LA][j]
	// I1 = DPI[LA][j]
		
		TBrow[j] = 0;
		float mi = Mrow[int(j)-1] + g_Viterbi_TermGapOpen;
		I1 += g_Viterbi_TermGapExt;
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
