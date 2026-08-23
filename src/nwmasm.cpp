#include "muscle.h"
#include "mx.h"
#include "xdpmem.h"
#include "masm.h"

void WriteLocalAln_MASM(FILE *f, const string &LabelA, const MASM &MA,
  const string &LabelQ, const vector<vector<byte> > &Q,
  uint Loi, uint Loj, const char *Path);

float NWFast_MASM_MegaProf(XDPMem &Mem, const MASM &MA,
  const vector<vector<byte> > &PB, uint &Loj, string &Path);

void cmd_nwmasm()
	{
	const string &MasmFN = g_Arg1;
	const string &StructsFN = opt(query);

	FILE *fOut = CreateStdioFile(opt(output));

	Mega::RejectLegacyMega(StructsFN);
	Mega::FromStructs(StructsFN);

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
		string Path;
		uint Loj;
		float Score = NWFast_MASM_MegaProf(Mem, M, Q, Loj, Path);
		WriteLocalAln_MASM(g_fLog, LabelM, M, LabelQ, Q, 0, Loj, Path.c_str());
		Log("Score = %.3g\n", Score);
		Log("\n");

		if (fOut != 0)
			{
			fprintf(fOut, "%s", LabelM.c_str());
			fprintf(fOut, "\t%s", LabelQ.c_str());
			fprintf(fOut, "\t%.3g", Score);
			fprintf(fOut, "\t%s", Path.c_str());
			fprintf(fOut, "\n");
			}
		}

	CloseStdioFile(fOut);
	}
