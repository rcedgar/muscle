#include "myutils.h"
#include "swtester.h"

static float GapOpen = -1;
static float GapExt = -0.4f;

/***
@SCOREDIFF Enum_MASM_Mega Simple_MASM_Mega
A: VDA|KMY|STN
B: MHS
  4.48/0.386  loa 5,1  lob 0,0  MMDM,MM
***/
static void Bug()
	{
	SWTester ST;

	SWer_Enum_MASM_Mega S;
	S.m_GapOpen = GapOpen;
	S.m_GapExt = GapExt;
	ST.RunXAB(S, "VDA|KMY|STN", "MHS", true);

	SWer_Simple_MASM_Mega S2;
	S2.m_GapOpen = GapOpen;
	S2.m_GapExt = GapExt;
	ST.RunXAB(S2, "VDA|KMY|STN", "MHS", true);
	}

static void TestSeqs()
	{
	SWTester ST;

	SWer_Simple_MASM_Mega Simple_MASM_Mega;
	SWer_MASM_Mega_Seqs MASM_Mega_Seqs;

	Simple_MASM_Mega.m_GapOpen = GapOpen;
	Simple_MASM_Mega.m_GapExt = GapExt;

	MASM_Mega_Seqs.m_GapOpen = GapOpen;
	MASM_Mega_Seqs.m_GapExt = GapExt;

	const uint MinL = 3;
	const uint MaxL = 5;
	const uint Iters = 1000;

	ST.ClearStats();
	ST.SetX(Simple_MASM_Mega);
	ST.SetY(MASM_Mega_Seqs);
	ST.RunRandomSeqsIters(MinL, MaxL, Iters);
	ST.Stats();
	}

void cmd_swtestmm()
	{
	//TestSeqs();
	//Bug();
	//return;

	SWTester ST;

	SWer_Enum_MASM_Mega Enum_MASM_Mega;
	SWer_Simple_MASM_Mega Simple_MASM_Mega;
	SWer_MASM_Mega MASM_Mega;

	Simple_MASM_Mega.m_GapOpen = GapOpen;
	Simple_MASM_Mega.m_GapExt = GapExt;

	MASM_Mega.m_GapOpen = GapOpen;
	MASM_Mega.m_GapExt = GapExt;

	Enum_MASM_Mega.m_GapOpen = GapOpen;
	Enum_MASM_Mega.m_GapExt = GapExt;

	const uint MinN = 1;
	const uint MaxN = 5;
	const uint MinL = 3;
	const uint MaxL = 7;
	const uint Iters = 1000;

	//ST.ClearStats();
	//ST.SetX(Simple_MASM_Mega);
	//ST.SetY(MASM_Mega);
	//ST.RunRandomMSASeqIters(MinN, MaxN, MinL, MaxL, Iters);
	//ST.Stats();

	ST.ClearStats();
	ST.SetX(Enum_MASM_Mega);
	ST.SetY(Simple_MASM_Mega);
	ST.RunRandomMSASeqIters(MinN, MaxN, MinL, MaxL, Iters);
	ST.Stats();
	}
