#include "myutils.h"
#include "swtester.h"

static float GapOpen = -1;
static float GapExt = -0.4f;

/***
@SCOREDIFF SWer_Simple_MASM_Mega Fast_Seqs_AA_BLOSUM62
A: EVRDYIQ
B: RQGEG
  3.98/3.58  loa 2,2  lob 0,0  MDDDM,MDDDM
***/
static void Bug()
	{
	SWer_Simple_MASM_Mega S;
	S.m_GapOpen = GapOpen;
	S.m_GapExt = GapExt;

	SWTester ST;
	ST.RunXAB(S, "EVRDYIQ", "RQGEG", true);

	SWer_Fast_Seqs_AA_BLOSUM62 S2;
	S2.m_GapOpen = GapOpen;
	S2.m_GapExt = GapExt;
	ST.RunXAB(S2, "EVRDYIQ", "RQGEG", true);
	}

void cmd_swtest()
	{
	SWTester ST;

	SWer_Enum_Seqs_AA_BLOSUM62 Enum_Seqs_AA_BLOSUM62;
	SWer_Fast_Seqs_AA_BLOSUM62 Fast_Seqs_AA_BLOSUM62;
	SWer_Simple_Seqs_AA_BLOSUM62 Simple_Seqs_AA_BLOSUM62;
	SWer_MASM_Mega_Seqs Mega_Prof_Seqs;
	SWer_Simple_MASM_Mega Simple_MASM_Mega;

	SWer_PS SPS; 
	PathScorer_AA_BLOSUM62 PSAB;
	PSAB.m_GapOpen = GapOpen;
	PSAB.m_GapExt = GapExt;
	SPS.m_PS = &PSAB;

	Enum_Seqs_AA_BLOSUM62.SetGaps(GapOpen, GapExt);

	Fast_Seqs_AA_BLOSUM62.m_GapOpen = GapOpen;
	Fast_Seqs_AA_BLOSUM62.m_GapExt = GapExt;

	Simple_Seqs_AA_BLOSUM62.m_GapOpen = GapOpen;
	Simple_Seqs_AA_BLOSUM62.m_GapExt = GapExt;

	Mega_Prof_Seqs.m_GapOpen = GapOpen;
	Mega_Prof_Seqs.m_GapExt = GapExt;

	Simple_MASM_Mega.m_GapOpen = GapOpen;
	Simple_MASM_Mega.m_GapExt = GapExt;

	const uint MinL = 3;
	const uint MaxL = 9;
	const uint Iters = 1000;
	//ST.RunXAB(Simple_MASM_Mega, "ESK", "MMMW", true);

	ST.ClearStats();
	ST.SetY(SPS);
	ST.SetX(Fast_Seqs_AA_BLOSUM62);
	ST.RunAB("SEQV", "EQ");
	//ST.RunRandomSeqsIters(MinL, MaxL, Iters);
	ST.Stats();
	return;

	ST.ClearStats();
	ST.SetY(Fast_Seqs_AA_BLOSUM62);
	ST.SetX(Simple_MASM_Mega);
	ST.RunRandomSeqsIters(MinL, MaxL, Iters);
	ST.Stats();

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

	ST.ClearStats();
	ST.SetY(Fast_Seqs_AA_BLOSUM62);
	ST.SetX(Mega_Prof_Seqs);
	ST.RunRandomSeqsIters(MinL, MaxL, Iters);
	ST.Stats();
	}
