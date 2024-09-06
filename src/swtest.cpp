#include "myutils.h"
#include "swtester.h"

void cmd_swtest()
	{
	SWTester ST;

	SWer_Enum_Seqs_AA_BLOSUM62 Enum_Seqs_AA_BLOSUM62;
	SWer_Fast_Seqs_AA_BLOSUM62 Fast_Seqs_AA_BLOSUM62;
	SWer_Simple_Seqs_AA_BLOSUM62 Simple_Seqs_AA_BLOSUM62;

	float GapOpen = -1;
	float GapExt = -0.5;

	Enum_Seqs_AA_BLOSUM62.SetGaps(GapOpen, GapExt);

	Fast_Seqs_AA_BLOSUM62.m_GapOpen = GapOpen;
	Fast_Seqs_AA_BLOSUM62.m_GapExt = GapExt;

	Simple_Seqs_AA_BLOSUM62.m_GapOpen = GapOpen;
	Simple_Seqs_AA_BLOSUM62.m_GapExt = GapExt;

	//ST.RunXAB(Enum_Seqs_AA_BLOSUM62, "SEQVE", "EQV");
	//ST.RunXY(Enum_Seqs_AA_BLOSUM62, Fast_Seqs_AA_BLOSUM62, "SEQVE", "EQV");
	//ST.RunXY(Enum_Seqs_AA_BLOSUM62, Simple_Seqs_AA_BLOSUM62, "SEQVE", "EQV");
	//ST.Stats();

	const uint MinL = 3;
	const uint MaxL = 7;
	const uint Iters = 1000;

	ST.ClearStats();
	ST.SetY(Fast_Seqs_AA_BLOSUM62);
	ST.SetX(Simple_Seqs_AA_BLOSUM62);
	ST.RunRandomSeqsIters(MinL, MaxL, Iters);
	ST.Stats();

	ST.ClearStats();
	ST.SetY(Fast_Seqs_AA_BLOSUM62);
	ST.SetX(Enum_Seqs_AA_BLOSUM62);
	ST.RunRandomSeqsIters(MinL, MaxL, Iters);
	ST.Stats();
	}
