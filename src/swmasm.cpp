#include "muscle.h"
#include "mx.h"
#include "tracebit.h"
#include "xdpmem.h"
#include "swtrace.h"
#include "masm.h"

void WriteLocalAln_MASM(FILE *f, const string &LabelA, const MASM &MA,
  const string &LabelQ, const vector<vector<byte> > &Q,
  uint Loi, uint Loj, const char *Path);

float SWFast_SMx(XDPMem &Mem, const Mx<float> &SMx,
  float Open, float Ext, uint &Loi, uint &Loj, uint &Leni, uint &Lenj,
  string &Path);

float SWFast_MASM(XDPMem &Mem, const MASM &A, const vector<vector<byte> > &B,
  uint &Loi, uint &Loj, uint &Leni, uint &Lenj, string &Path)
	{
	float Open = A.m_GapOpen;
	float Ext = A.m_GapExt;
	Mx<float> SMx;
	A.MakeSMx(B, SMx);
	float Score = SWFast_SMx(Mem, SMx, -Open, -Ext, Loi, Loj, Leni, Lenj, Path);
	return Score;
	}

void cmd_swmasm()
	{
	const string &MasmFN = g_Arg1;
	const string &MegaFN = opt(query);

	FILE *fOut = CreateStdioFile(opt(output));

	Mega::FromFile(MegaFN);

	MASM M;
	M.FromFile(MasmFN);
	const string &LabelM = M.m_Label;

	XDPMem Mem;
	const uint QueryProfileCount = Mega::GetProfileCount();
	for (uint i = 0; i < QueryProfileCount; ++i)
		{
		ProgressStep(i, QueryProfileCount, "Aligning");
		const vector<vector<byte> > &Q = Mega::GetProfile(i);
		const string &LabelQ = Mega::GetLabel(i);
		uint Loi, Loj, Leni, Lenj;
		string Path;
		float Score = SWFast_MASM(Mem, M, Q, Loi, Loj, Leni, Lenj, Path);
		WriteLocalAln_MASM(g_fLog, LabelM, M, LabelQ, Q, Loi, Loj, Path.c_str());
		Log("Score = %.3g\n", Score);
		Log("\n");

		if (fOut != 0)
			{
			fprintf(fOut, "%s", LabelM.c_str());
			fprintf(fOut, "\t%s", LabelQ.c_str());
			fprintf(fOut, "\t%.3g", Score);
			fprintf(fOut, "\n");
			}
		}

	CloseStdioFile(fOut);
	}
