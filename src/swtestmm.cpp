#include "myutils.h"
#include "swtester.h"

static float GapOpen = -1;
static float GapExt = -0.4f;

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
	TestSeqs();

	SWTester ST;

	SWer_Simple_MASM_Mega Simple_MASM_Mega;
	SWer_MASM_Mega MASM_Mega;

	Simple_MASM_Mega.m_GapOpen = GapOpen;
	Simple_MASM_Mega.m_GapExt = GapExt;

	MASM_Mega.m_GapOpen = GapOpen;
	MASM_Mega.m_GapExt = GapExt;

	const uint MinN = 1;
	const uint MaxN = 5;
	const uint MinL = 3;
	const uint MaxL = 5;
	const uint Iters = 1000;

	ST.ClearStats();
	ST.SetX(Simple_MASM_Mega);
	ST.SetY(MASM_Mega);
	ST.RunRandomMSASeqIters(MinN, MaxN, MinL, MaxL, Iters);
	ST.Stats();
	}
