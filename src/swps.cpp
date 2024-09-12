#include "muscle.h"
#include "mx.h"
#include "tracebit.h"
#include "xdpmem.h"
#include "pathscorer.h"
#include "allocmx.h"

void TraceBackBitSW(XDPMem &Mem,
  uint LA, uint LB, uint Besti, uint Bestj,
  uint &Leni, uint &Lenj, string &Path);

float SWPS(XDPMem &Mem, PathScorer &PS, uint &Loi, uint &Loj, string &Path)
	{
	const uint LA = PS.GetLA();
	const uint LB = PS.GetLB();
	Mem.Alloc(LA+32, LB+32);

	float *Mrow = Mem.GetDPRow1();
	float *Drow = Mem.GetDPRow2();
	byte **TB = Mem.GetTBBit();

	for (uint j = 0; j <= LB; ++j)
		{
		Mrow[j] = MINUS_INFINITY;
		Drow[j] = MINUS_INFINITY;
		}
	
	float BestScore = 0;
	uint Besti = UINT_MAX;
	uint Bestj = UINT_MAX;

	for (uint i = 0; i < LA; ++i)
		{
		float Is = MINUS_INFINITY;
		byte *TBrow = TB[i];
		float PrevSavedM = MINUS_INFINITY;
		for (uint j = 0; j < LB; ++j)
			{
			byte TraceBits = 0;
			float SavedM = Mrow[j+1];

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

			float mm = PrevSavedM + PS.GetScoreMM(i, j);
			float dm = Drow[j] + PS.GetScoreDM(i, j);
			float im = Is + PS.GetScoreIM(i, j);

			////float t = max(max(mm, dm), im);
			////if (t < 0)
			////	t = 0;
			////float m = t + PS.GetMatchScore(i, j);
			////BestScore = max(m, BestScore);

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
				t = 0;
				TraceBits = TRACEBITS_SM;
				}
			float m = t + PS.GetMatchScore(i, j);
			if (m > BestScore)
				{
				BestScore = m;
				Besti = i;
				Bestj = j;
				}

			// Match output:
			//		Mrow[j+1] = M[I][j+1]
			Mrow[j+1] = m;

			////////////////////////////
			//	D[i+1][j] = max
			//		M[i][j] + tMD,
			//		D[i][j] + tDD
			////////////////////////////

			// Delete input:
			//		PrevSavedM = M[i][j]
			//		Drow[j] = D[i][j]
			float md = PrevSavedM + PS.GetScoreMD(i, j);
			float dd = Drow[j] + PS.GetScoreDD(i, j);

			float d = md;
			if (dd > d)
				{
				d = md;
				TraceBits |= TRACEBITS_MD;
				}

			// Delete output:
			//		Drow[j] = D[I][j]
			Drow[j] = d;

			// Insert input:
			//	PrevSavedM = M[i][j]
			//	Is = I[i][j]
			float mi = PrevSavedM + PS.GetScoreMI(i, j);
			float ii = Is + PS.GetScoreII(i, j);

			// Insert output:
			//	Is = I[i][j+1]
			////Is = max(mi, ii);
			Is = ii;
			if (mi > Is)
				{
				Is = mi;
				TraceBits |= TRACEBITS_MI;
				}

			PrevSavedM = SavedM;
			}
		}

	if (BestScore <= 0.0f)
		return 0.0f;

	uint Leni, Lenj;
	TraceBackBitSW(Mem, LA, LB, Besti+1, Bestj+1,
	  Leni, Lenj, Path);
	asserta(Besti+1 >= Leni);
	asserta(Bestj+1 >= Lenj);

	Loi = Besti + 1 - Leni;
	Loj = Bestj + 1 - Lenj;

	return BestScore;
	}
