#include "muscle.h"
#include "mx.h"
#include "tracebit.h"
#include "xdpmem.h"
#include "swtrace.h"
#include "masm.h"

void WriteLocalAln(FILE *f, const string &LabelA, const byte *A,
  const string &LabelB, const byte *B,
  uint Loi, uint Loj, const char *Path);
float SWFast_SMx(XDPMem &Mem, const Mx<float> &SMx,
  float Open, float Ext, uint &Loi, uint &Loj, uint &Leni, uint &Lenj,
  string &Path);

float SWFast_MASM_Seq(XDPMem &Mem, const MASM &A, const Sequence &B,
  float Open, float Ext,
  uint &Loi, uint &Loj, uint &Leni, uint &Lenj, string &Path)
	{
	Mx<float> SMx;
	A.MakeSMx_Sequence(B, SMx);
	float Score = SWFast_SMx(Mem, SMx, -Open, -Ext, Loi, Loj, Leni, Lenj, Path);
	return Score;
	}

void cmd_swmasm_seq()
	{
	const string &AlnFN = g_Arg1;
	const string &MegaFN = opt(input);
	const string &FaFN = opt(input2);

	Mega::FromFile(MegaFN);

	MultiSequence Aln;
	Aln.FromFASTA(AlnFN);

	float GapOpen = 4;
	float GapExt = 0.5;

	MASM M;
	M.FromMSA(Aln, "FomMSA", GapOpen, GapExt);
	M.ToFile(opt(output));

	MultiSequence Query;
	Query.FromFASTA(FaFN);

	XDPMem Mem;
	const uint QuerySeqCount = Query.GetSeqCount();
	for (uint i = 0; i < QuerySeqCount; ++i)
		{
		const Sequence &Q = *Query.GetSequence(i);
		uint Loi, Loj, Leni, Lenj;
		string Path;
		float Score = SWFast_MASM_Seq(Mem, M, Q, GapOpen, GapExt,
		  Loi, Loj, Leni, Lenj, Path);
		Log("%10.3g  %16.16s  %7u  %7u  %s\n",
		  Score, Q.GetLabel().c_str(), Loi, Loj, Path.c_str());
		Log("\n");
		}
	}
