#include "muscle.h"
#include "mx.h"
#include "xdpmem.h"
#include "masm.h"
#include "sequence.h"

void WriteLocalAln(FILE *f, const string &LabelA, const byte *A,
  const string &LabelB, const byte *B,
  uint Loi, uint Loj, const char *Path);

float NWFast_MASM_Seq(XDPMem &Mem, const MASM &A, const Sequence &B,
  string &Path);

void cmd_nwmasm_seq()
	{
	const string &AlnFN = g_Arg1;
	const string &StructsFN = opt(input);
	const string &FaFN = opt(input2);

	Mega::RejectLegacyMega(StructsFN);
	Mega::FromStructs(StructsFN);

	MultiSequence Aln;
	Aln.FromFASTA(AlnFN);

	float GapOpen = 4;
	float GapExt = 0.5;
	if (optset_gapopen)
		GapOpen = (float) opt(gapopen);
	if (optset_gapext)
		GapExt = (float) opt(gapext);

	MASM M;
	M.FromMSA(Aln, "FromMSA", GapOpen, GapExt);
	M.ToFile(opt(output));

	MultiSequence Query;
	Query.FromFASTA(FaFN);

	XDPMem Mem;
	const uint QuerySeqCount = Query.GetSeqCount();
	string Cons;
	M.GetConsensusSeq(Cons);
	for (uint i = 0; i < QuerySeqCount; ++i)
		{
		const Sequence &Q = *Query.GetSequence(i);
		string Path;
		float Score = NWFast_MASM_Seq(Mem, M, Q, Path);
		WriteLocalAln(g_fLog, M.m_Label.c_str(), (const byte *) Cons.c_str(),
		  Q.GetLabelCStr(), Q.GetBytePtr(),
		  0, 0, Path.c_str());
		Log("%10.3g  %16.16s  %s\n",
		  Score, Q.GetLabel().c_str(), Path.c_str());
		Log("\n");
		}
	}
