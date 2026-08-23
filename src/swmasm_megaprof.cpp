#include "muscle.h"
#include "mx.h"
#include "tracebit.h"
#include "xdpmem.h"
#include "swtrace.h"
#include "masm.h"
#include "sequence.h"
#include <algorithm>

void TraceBackBitSW(XDPMem &Mem,
  uint LA, uint LB, uint Besti, uint Bestj,
  uint &Leni, uint &Lenj, string &Path);

static float InsertOpen(const MASM &MA)
	{
	asserta(MA.m_GapOpen != FLT_MAX);
	return -MA.m_GapOpen/2;
	}

static float InsertExt(const MASM &MA)
	{
	asserta(MA.m_GapExt != FLT_MAX);
	return -MA.m_GapExt;
	}

static void AssertColGaps(const MASMCol &Col)
	{
	asserta(Col.m_GapOpen != FLT_MAX);
	asserta(Col.m_GapExt != FLT_MAX);
	asserta(Col.m_GapClose != FLT_MAX);
	}

static void TraceBackBitNW(XDPMem &Mem, uint LA, uint Startj, char State,
  uint &Loj, string &Path)
	{
	Path.clear();
	byte **TB = Mem.GetTBBit();
	uint i = LA;
	uint j = Startj;
	for (;;)
		{
		if (i == 0)
			break;
		Path += State;
		byte t;
		switch (State)
			{
		case 'M':
			asserta(i > 0 && j > 0);
			t = TB[i-1][j-1];
			if (t & TRACEBITS_DM)
				State = 'D';
			else if (t & TRACEBITS_IM)
				State = 'I';
			else
				State = 'M';
			--i;
			--j;
			break;
		case 'D':
			asserta(i > 0);
			t = TB[i-1][j];
			if (t & TRACEBITS_MD)
				State = 'M';
			else
				State = 'D';
			--i;
			break;
		case 'I':
			asserta(j > 0);
			t = TB[i][j-1];
			if (t & TRACEBITS_MI)
				State = 'M';
			else
				State = 'I';
			--j;
			break;
		default:
			Die("TraceBackBitNW, invalid state %c", State);
			}
		}
	Loj = j;
	reverse(Path.begin(), Path.end());
	}

static float SWFast_MASM_SMx(XDPMem &Mem, const MASM &MA, const Mx<float> &SMx,
  uint &Loi, uint &Loj, uint &Leni, uint &Lenj, string &Path)
	{
	const uint LA = MA.GetColCount();
	const uint LB = SMx.GetColCount();
	asserta(SMx.GetRowCount() == LA);
	asserta(LA > 0 && LB > 0);

	const float OpenI = InsertOpen(MA);
	const float ExtI = InsertExt(MA);
	const float * const *SMxData = SMx.GetData();

	Mem.Alloc(LA+32, LB+32);

	Leni = 0;
	Lenj = 0;

	float *Mrow = Mem.GetDPRow1();
	float *Drow = Mem.GetDPRow2();
	byte **TB = Mem.GetTBBit();
	INIT_TRACE(LA, LB, TB);

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

	float M0 = 0.0f;
	for (uint i = 0; i < LA; ++i)
		{
		const MASMCol &ColA = MA.GetCol(i);
		AssertColGaps(ColA);
		const float OpenD = -ColA.m_GapOpen;
		const float ExtD = -ColA.m_GapExt;
		const float CloseD = -ColA.m_GapClose;
		const float *SMxRow = SMxData[i];
		float I0 = MINUS_INFINITY;
		byte *TBrow = TB[i];
		for (uint j = 0; j < LB; ++j)
			{
			byte TraceBits = 0;
			float SavedM0 = M0;

			float xM = M0;
			if (Drow[j] + CloseD > xM)
				{
				xM = Drow[j] + CloseD;
				TraceBits = TRACEBITS_DM;
				}
			if (I0 + OpenI > xM)
				{
				xM = I0 + OpenI;
				TraceBits = TRACEBITS_IM;
				}
			if (0.0f >= xM)
				{
				xM = 0.0f;
				TraceBits = TRACEBITS_SM;
				}

			M0 = Mrow[j];
			xM += SMxRow[j];
			if (xM > BestScore)
				{
				BestScore = xM;
				Besti = i;
				Bestj = j;
				}

			Mrow[j] = xM;
			TRACE_M(i, j, xM);

			float md = SavedM0 + OpenD;
			Drow[j] += ExtD;
			if (md >= Drow[j])
				{
				Drow[j] = md;
				TraceBits |= TRACEBITS_MD;
				}
			TRACE_D(i, j, Drow[j]);

			float mi = SavedM0 + OpenI;
			I0 += ExtI;
			if (mi >= I0)
				{
				I0 = mi;
				TraceBits |= TRACEBITS_MI;
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

static float NWFast_MASM_SMx(XDPMem &Mem, const MASM &MA, const Mx<float> &SMx,
  uint &Loj, string &Path)
	{
	const uint LA = MA.GetColCount();
	const uint LB = SMx.GetColCount();
	asserta(SMx.GetRowCount() == LA);
	asserta(LA > 0 && LB > 0);

	const float OpenI = InsertOpen(MA);
	const float ExtI = InsertExt(MA);
	const float * const *SMxData = SMx.GetData();

	Mem.Alloc(LA+32, LB+32);

	float *Mrow = Mem.GetDPRow1();
	float *Drow = Mem.GetDPRow2();
	byte **TB = Mem.GetTBBit();

	Mrow[-1] = MINUS_INFINITY;
	for (uint j = 0; j <= LB; ++j)
		{
		Mrow[j] = MINUS_INFINITY;
		Drow[j] = MINUS_INFINITY;
		}

	float M0 = 0.0f;
	for (uint i = 0; i < LA; ++i)
		{
		const MASMCol &ColA = MA.GetCol(i);
		AssertColGaps(ColA);
		const float OpenD = -ColA.m_GapOpen;
		const float ExtD = -ColA.m_GapExt;
		const float CloseD = -ColA.m_GapClose;
		const float *SMxRow = SMxData[i];
		float I0 = MINUS_INFINITY;
		byte *TBrow = TB[i];
		for (uint j = 0; j < LB; ++j)
			{
			byte TraceBits = 0;
			float SavedM0 = M0;

			float xM;
			if (i == 0)
				{
			// Free leading query gaps: match model col 0 at any j
				xM = 0.0f;
				}
			else
				{
				xM = M0;
				if (Drow[j] + CloseD > xM)
					{
					xM = Drow[j] + CloseD;
					TraceBits = TRACEBITS_DM;
					}
				if (I0 + OpenI > xM)
					{
					xM = I0 + OpenI;
					TraceBits = TRACEBITS_IM;
					}
				}

			M0 = Mrow[j];
			Mrow[j] = xM + SMxRow[j];

			if (i == 0)
				{
				Drow[j] = OpenD;
				TraceBits |= TRACEBITS_MD;
				}
			else
				{
				float md = SavedM0 + OpenD;
				Drow[j] += ExtD;
				if (md >= Drow[j])
					{
					Drow[j] = md;
					TraceBits |= TRACEBITS_MD;
					}
				}

			float mi = SavedM0 + OpenI;
			I0 += ExtI;
			if (mi >= I0)
				{
				I0 = mi;
				TraceBits |= TRACEBITS_MI;
				}

			TBrow[j] = TraceBits;
			}

		TBrow[LB] = 0;
		float md = M0 + OpenD;
		Drow[LB] += ExtD;
		if (md >= Drow[LB])
			{
			Drow[LB] = md;
			TBrow[LB] = TRACEBITS_MD;
			}

		M0 = MINUS_INFINITY;
		}

	float Score = MINUS_INFINITY;
	char State = 'M';
	uint Endj = 0;
	for (uint j = 0; j < LB; ++j)
		{
		if (Mrow[j] > Score)
			{
			Score = Mrow[j];
			State = 'M';
			Endj = j;
			}
		if (Drow[j] > Score)
			{
			Score = Drow[j];
			State = 'D';
			Endj = j;
			}
		}
	if (Drow[LB] > Score)
		{
		Score = Drow[LB];
		State = 'D';
		Endj = LB;
		}

	const uint Startj = (State == 'M' ? Endj + 1 : Endj);
	TraceBackBitNW(Mem, LA, Startj, State, Loj, Path);
	return Score;
	}

float SWFast_MASM_MegaProf(XDPMem &Mem, const MASM &MA,
  const vector<vector<byte> > &PB,
  uint &Loi, uint &Loj, uint &Leni, uint &Lenj, string &Path)
	{
	Mx<float> SMx;
	MA.MakeSMx(PB, SMx);
	return SWFast_MASM_SMx(Mem, MA, SMx, Loi, Loj, Leni, Lenj, Path);
	}

float NWFast_MASM_MegaProf(XDPMem &Mem, const MASM &MA,
  const vector<vector<byte> > &PB, uint &Loj, string &Path)
	{
	Mx<float> SMx;
	MA.MakeSMx(PB, SMx);
	return NWFast_MASM_SMx(Mem, MA, SMx, Loj, Path);
	}

float SWFast_MASM(XDPMem &Mem, const MASM &A, const vector<vector<byte> > &B,
  uint &Loi, uint &Loj, uint &Leni, uint &Lenj, string &Path)
	{
	return SWFast_MASM_MegaProf(Mem, A, B, Loi, Loj, Leni, Lenj, Path);
	}

float SWFast_MASM_Seq(XDPMem &Mem, const MASM &A, const Sequence &B,
  uint &Loi, uint &Loj, uint &Leni, uint &Lenj, string &Path)
	{
	Mx<float> SMx;
	A.MakeSMx_Sequence(B, SMx);
	return SWFast_MASM_SMx(Mem, A, SMx, Loi, Loj, Leni, Lenj, Path);
	}

float NWFast_MASM_Seq(XDPMem &Mem, const MASM &A, const Sequence &B,
  uint &Loj, string &Path)
	{
	Mx<float> SMx;
	A.MakeSMx_Sequence(B, SMx);
	return NWFast_MASM_SMx(Mem, A, SMx, Loj, Path);
	}
