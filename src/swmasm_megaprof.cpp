#include "muscle.h"
#include "mx.h"
#include "tracebit.h"
#include "xdpmem.h"
#include "swtrace.h"
#include "masm.h"

void TraceBackBitSW(XDPMem &Mem,
  uint LA, uint LB, uint Besti, uint Bestj,
  uint &Leni, uint &Lenj, string &Path);

float SWFast_MASM_MegaProf(XDPMem &Mem, const MASM &MA,
  const vector<vector<byte> > &PB, float Open, float Ext,
  uint &Loi, uint &Loj, uint &Leni, uint &Lenj, string &Path)
	{
#if TRACE && !DOTONLY
	SMx.LogMe();
#endif
	const uint LA = MA.GetColCount();
	const uint LB = SIZE(PB);
	asserta(Open <= 0);
	asserta(Ext <= 0);

	Mem.Alloc(LA+32, LB+32);

	Leni = 0;
	Lenj = 0;

	float *Mrow = Mem.GetDPRow1();
	float *Drow = Mem.GetDPRow2();
	byte **TB = Mem.GetTBBit();
	INIT_TRACE(LA, LB, TB);

// Use Mrow[-1], so...
	Mrow[-1] = MINUS_INFINITY;
	TRACE_M(0, -1, MINUS_INFINITY);

	for (uint j = 0; j <= LB; ++j)
		{
		Mrow[j] = MINUS_INFINITY;
		Drow[j] = MINUS_INFINITY;
		TRACE_M(0, j, MINUS_INFINITY);
		TRACE_D(0, j, MINUS_INFINITY);
		}
	
	float BestScore = 0.0f;
	uint Besti = UINT_MAX;
	uint Bestj = UINT_MAX;

// Main loop
	float M0 = float (0);
	for (uint i = 0; i < LA; ++i)
		{
		//const float *SMxRow = SMxData[i];
		const MASMCol &ColA = MA.GetCol(i);
		float I0 = MINUS_INFINITY;
		byte *TBrow = TB[i];
		for (uint j = 0; j < LB; ++j)
			{
			const vector<byte> &ColB = PB[j];
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
			if (0.0f >= xM)
				{
				xM = 0.0f;
				TraceBits = TRACEBITS_SM;
				}

			M0 = Mrow[j];
			float MatchScore = ColA.GetMatchScore_MegaProfilePos(ColB);
			//xM += SMxRow[j];
			xM += MatchScore;
			if (xM > BestScore)
				{
				BestScore = xM;
				Besti = i;
				Bestj = j;
				}

			Mrow[j] = xM;
			TRACE_M(i, j, xM);
		// Mrow[j] = DPM[i+1][j+1])
			}
			
		// DELETE
			{
		// SavedM0 = DPM[i][j]
		// Drow[j] = DPD[i][j]
			float md = SavedM0 + Open;
			Drow[j] += Ext;
			if (md >= Drow[j])
				{
				Drow[j] = md;
				TraceBits |= TRACEBITS_MD;
				}
			TRACE_D(i, j, Drow[j]);
		// Drow[j] = DPD[i+1][j]
			}
			
		// INSERT
			{
		// SavedM0 = DPM[i][j]
		// I0 = DPI[i][j]
			float mi = SavedM0 + Open;
			I0 += Ext;
			if (mi >= I0)
				{
				I0 = mi;
				TraceBits |= TRACEBITS_MI;
				}
			}
			
			TBrow[j] = TraceBits;
			}
		
		M0 = MINUS_INFINITY;
		}

	DONE_TRACE(BestScore, Besti, Bestj, TB);
	if (BestScore <= 0.0f)
		return 0.0f;

	TraceBackBitSW(Mem, LA, LB, Besti+1, Bestj+1,
	  Leni, Lenj, Path);
	asserta(Besti+1 >= Leni);
	asserta(Bestj+1 >= Lenj);

	Loi = Besti + 1 - Leni;
	Loj = Bestj + 1 - Lenj;

	return BestScore;
	}
