#include "muscle.h"
#include "mx.h"
#include "tracebit.h"
#include "xdpmem.h"
#include "swtrace.h"
#include "masm.h"

void WriteLocalAln_MASM(FILE *f, const MASM &MA, const vector<vector<byte> > &Q,
  uint Loi, uint Loj, const char *Path);

float SWFast_SMx(XDPMem &Mem, const Mx<float> &SMx,
  float Open, float Ext, uint &Loi, uint &Loj, uint &Leni, uint &Lenj,
  string &Path);

float SWFast_MASM(XDPMem &Mem, const MASM &A, const vector<vector<byte> > &B,
  float Open, float Ext,
  uint &Loi, uint &Loj, uint &Leni, uint &Lenj, string &Path)
	{
	Mx<float> SMx;
	A.MakeSMx(B, SMx);
	float Score = SWFast_SMx(Mem, SMx, -Open, -Ext, Loi, Loj, Leni, Lenj, Path);
	return Score;
	}

void cmd_swmasm()
	{
	const string &MasmFN = g_Arg1;
	const string &MegaFN = opt(query);

	Mega::FromFile(MegaFN);

	MASM M;
	M.FromFile(MasmFN);

	float GapOpen = 1.5;
	float GapExt = 0.42;

	XDPMem Mem;
	const uint QueryProfileCount = Mega::GetProfileCount();
	for (uint i = 0; i < QueryProfileCount; ++i)
		{
		const vector<vector<byte> > &Q = Mega::GetProfile(i);
		const string &LabelQ = Mega::GetLabel(i);
		uint Loi, Loj, Leni, Lenj;
		string Path;
		float Score = SWFast_MASM(Mem, M, Q, GapOpen, GapExt,
		  Loi, Loj, Leni, Lenj, Path);
		WriteLocalAln_MASM(g_fLog, M, Q, Loi, Loj, Path.c_str());
		Log("%10.3g  %16.16s  %7u  %7u  %s\n",
		  Score, LabelQ.c_str(), Loi, Loj, Path.c_str());
		Log("\n");
		}
	}
