#include "muscle.h"
#include "mx.h"
#include "tracebit.h"
#include "xdpmem.h"
#include "swtrace.h"

void MakeBlosum62SMx(const Sequence &A, const Sequence &B, Mx<float> &MxS);
void MakeBlosum62SMx(const string &A, const string &B, Mx<float> &MxS);

void LogTBSW(const char *Msg, XDPMem &Mem, uint LA, uint LB)
	{
	byte **TB = Mem.GetTBBit();
	Log("TBM %s\n", Msg);
	for (uint i = 0; i < LA; ++i)
		{
		Log("%3u | ", i);
		for (uint j = 0; j < LB; ++j)
			{
			char c = 'M';
			byte Bits = TB[i][j];
			if (Bits == TRACEBITS_UNINIT)
				c = '*';
			else if (Bits & TRACEBITS_DM)
				c = 'D';
			else if (Bits & TRACEBITS_IM)
				c = 'I';
			else if (Bits & TRACEBITS_SM)
				c = 'S';
			Log("%c", c);
			}
		Log("\n");
		}

	Log("\nTBD %s\n", Msg);
	for (uint i = 0; i < LA; ++i)
		{
		Log("%3u | ", i);
		for (uint j = 0; j < LB; ++j)
			{
			char c = 'D';
			byte Bits = TB[i][j];
			if (Bits == TRACEBITS_UNINIT)
				c = '*';
			else if (Bits & TRACEBITS_MD)
				c = 'M';
			Log("%c", c);
			}
		Log("\n");
		}

	Log("\nTBI\n", Msg);
	for (uint i = 0; i < LA; ++i)
		{
		Log("%3u | ", i);
		for (uint j = 0; j < LB; ++j)
			{
			char c = 'I';
			byte Bits = TB[i][j];
			if (Bits == TRACEBITS_UNINIT)
				c = '*';
			else if (Bits & TRACEBITS_MI)
				c = 'M';
			Log("%c", c);
			}
		Log("\n");
		}
	}

void TraceBackBitSW(XDPMem &Mem,
  uint LA, uint LB, uint Besti, uint Bestj,
  uint &Leni, uint &Lenj, string &Path)
	{
	asserta(Besti < LA);
	asserta(Bestj < LB);
	Path.clear();
	byte **TB = Mem.GetTBBit();

#if TRACE && !DOTONLY
	Log("\n");
	Log("TraceBackBitSW\n");
#endif

	uint i = Besti;
	uint j = Bestj;
	char State = 'M';
	for (;;)
		{
#if TRACE && !DOTONLY
		Log("i=%3d  j=%3d  state=%c\n", (int) i, (int) j, State);
#endif
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
			else if (t & TRACEBITS_SM)
				{
				Leni = Besti - i + 1;
				Lenj = Bestj - j + 1;
				reverse(Path.begin(), Path.end());
				return;
				}
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
			Die("TraceBackBitSW, invalid state %c", State);
			}
		}
	}

float SWFast_SMx(XDPMem &Mem, const Mx<float> &SMx,
  float Open, float Ext, uint &Loi, uint &Loj, uint &Leni, uint &Lenj,
  string &Path)
	{
#if TRACE && !DOTONLY
	SMx.LogMe();
#endif
	const uint LA = SMx.GetRowCount();
	const uint LB = SMx.GetColCount();
	asserta(Open <= 0);
	asserta(Ext <= 0);

	Mem.Alloc(LA+32, LB+32);
	const float * const *SMxData = SMx.GetData();

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
	float M0 = 0;
	for (uint i = 0; i < LA; ++i)
		{
		const float *SMxRow = SMxData[i];
		float I0 = MINUS_INFINITY;
		byte *TBrow = TB[i];
		for (uint j = 0; j < LB; ++j)
			{
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
			xM += SMxRow[j];
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

float SWFast_Strings_BLOSUM62(XDPMem &Mem,
  const string &A, const string &B, float Open, float Ext,
  uint &Loi, uint &Loj, uint &Leni, uint &Lenj, string &Path)
	{
	Mx<float> SMx;
	MakeBlosum62SMx(A, B, SMx);
	float Score = SWFast_SMx(Mem, SMx, Open, Ext, Loi, Loj, Leni, Lenj, Path);
	return Score;
	}

float SWFast_Seqs_BLOSUM62(XDPMem &Mem,
  const Sequence &A, const Sequence &B, float Open, float Ext,
  uint &Loi, uint &Loj, uint &Leni, uint &Lenj, string &Path)
	{
	Mx<float> SMx;
	MakeBlosum62SMx(A, B, SMx);
	float Score = SWFast_SMx(Mem, SMx, Open, Ext, Loi, Loj, Leni, Lenj, Path);
	return Score;
	}

static void MakeLogOddsSMx(const Sequence &A, const Sequence &B, 
  const vector<vector<float> > &LogOddsMx, Mx<float> &MxS)
	{
	uint LA = A.GetLength();
	uint LB = B.GetLength();
	const char *ptrA = A.GetCharPtr();
	const char *ptrB = B.GetCharPtr();
	MxS.Alloc("LOS", LA, LB);
	float **S = MxS.GetData();
	for (uint i = 0; i < LA; ++i)
		{
		char a = ptrA[i];
		byte ai = g_CharToLetterAmino[a];
		for (uint j = 0; j < LB; ++j)
			{
			char b = ptrB[j];
			byte bi = g_CharToLetterAmino[b];
			if (ai < 20 && bi < 20)
				S[i][j] = LogOddsMx[ai][bi];
			else
				S[i][j] = 0;
			}
		}
	}

float SWFast_Seqs_LO(XDPMem &Mem, const vector<vector<float> > &LogOddsMx,
  const Sequence &A, const Sequence &B, float Open, float Ext,
  uint &Loi, uint &Loj, uint &Leni, uint &Lenj, string &Path)
	{
	Mx<float> SMx;
	MakeLogOddsSMx(A, B, LogOddsMx, SMx);
	float Score = SWFast_SMx(Mem, SMx, Open, Ext, Loi, Loj, Leni, Lenj, Path);
	return Score;
	}

void cmd_sw()
	{
	const string &FaFN = g_Arg1;

	MultiSequence Input;
	LoadInput(Input);

	const uint SeqCount = Input.GetSeqCount();
	if (SeqCount < 2)
		Die("%u seqs", SeqCount);

	float Open = -5;
	float Ext = -1;

	XDPMem Mem;
	uint PairCount = (SeqCount*(SeqCount - 1))/2;
	uint PairIdx = 0;
	for (uint i = 0; i < SeqCount; ++i)
		{
		const Sequence &seq_i = *Input.GetSequence(i);
		const string &Label_i = seq_i.m_Label;
		for (uint j = i+1; j < SeqCount; ++j)
			{
			ProgressStep(PairIdx++, PairCount, "Aligning");
			const Sequence &seq_j = *Input.GetSequence(j);
			const string &Label_j = seq_j.m_Label;

			uint Loi, Loj, Leni, Lenj;
			string Path;
			float Score = SWFast_Seqs_BLOSUM62(Mem, seq_i, seq_j, Open, Ext,
			  Loi, Loj, Leni, Lenj, Path);

			Log("\n");
			WriteLocalAln(g_fLog, Label_i, seq_i.GetBytePtr(),
			  Label_j, seq_j.GetBytePtr(), Loi, Loj, Path.c_str());
			Log("%s %s %.3g\n",
			  seq_i.GetLabel().c_str(),
			  seq_j.GetLabel().c_str(),
			  Score);
			}
		}
	}
