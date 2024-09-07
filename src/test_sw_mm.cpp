#include "muscle.h"
#include "pathscorer.h"
#include "enumpaths.h"
#include "xdpmem.h"

float SWFast_MASM_MegaProf(XDPMem &Mem, const MASM &MA,
  const vector<vector<byte> > &PB, float Open, float Ext,
  uint &Loi, uint &Loj, uint &Leni, uint &Lenj, string &Path);

void MakeMegaProfile_AA(const string &Seq,
  vector<vector<byte> > &Prof);

static uint g_MinRandL = 3;
static uint g_MaxRandL = 7;
static uint g_MinRandN = 1;
static uint g_MaxRandN = 5;
static uint g_RandGapPct = 25;

static uint g_N = 0;
static uint g_NAllOk = 0;
static uint g_NPathDiff = 0;
static uint g_NScoreDiff = 0;

static uint g_LA;
static uint g_LB;
static uint g_SeqCountA;
static float g_GapOpen = -3;
static float g_GapExt = -1;
static PathScorer_MASM_Mega *g_PS;
static float g_BestScore;
static string g_BestPath;
static uint g_BestPosA;
static uint g_BestPosB;
static bool g_LogAllPaths = false;

static void ClearBrute()
	{
	g_BestScore = 0;
	g_BestPath.clear();
	g_BestPosA = UINT_MAX;
	g_BestPosB = UINT_MAX;
	}

static void OnPath(uint PosA, uint PosB, const string &Path)
	{
	g_PS->m_LA = g_LA;
	g_PS->m_LB = g_LB;
	float Score = g_PS->GetLocalScore(PosA, PosB, Path);
	if (g_LogAllPaths)
		Log("%10.3g  %5u  %5u  %s\n", Score, PosA, PosB, Path.c_str());
	if (Score > g_BestScore)
		{
		g_BestScore = Score;
		g_BestPath = Path;
		g_BestPosA = PosA;
		g_BestPosB = PosB;
		}
	}

static MASM *MakeMASM_AAs(const vector<string> &Rows)
	{
	MultiSequence *Aln = new MultiSequence;
	const uint SeqCount = SIZE(Rows);
	uint ColCount = UINT_MAX;
	for (uint i = 0; i < SeqCount; ++i)
		{
		const string &Row = Rows[i];
		if (i == 0)
			ColCount = SIZE(Row);
		else
			asserta(SIZE(Row) == ColCount);
		Sequence *s = NewSequence();
		s->FromString("Row", Row);
		Aln->AddSequence(s, true);
		}
	MASM *M = new MASM;
	Mega::FromMSA_AAOnly(*Aln, g_GapOpen, g_GapExt);
	M->FromMSA(*Aln, "MSA", -g_GapOpen, -g_GapExt);
	return M;
	}

//static void Test_MASM_Mega(const vector<string> &RowsA, const string &B)
//	{
//	MASM &MA = *MakeMASM_AAs(RowsA);
//	vector<vector<byte> > PB;
//	MakeMegaProfile_AA(B, PB);
//	XDPMem Mem;
//	uint Loi, Loj, Leni, Lenj;
//	string Path;
//	float Score = SWFast_MASM_MegaProf(Mem, MA, PB, g_GapOpen, g_GapExt,
//	  Loi, Loj, Leni, Lenj, Path);
//	Log("Test_MASM_Mega %.3g (%u, %u) %s\n",
//	  Score, Loi, Loj, Path.c_str());
//	}

static void Test(const vector<string> &RowsA, const string &B,
  bool AlwaysLog)
	{
	g_SeqCountA = SIZE(RowsA);
	for (uint i = 0; i < g_SeqCountA; ++i)
		{
		if (i == 0)
			g_LA = SIZE(RowsA[i]);
		else
			asserta(SIZE(RowsA[i]) == g_LA);
		}

	XDPMem Mem;
	MASM &MA = *MakeMASM_AAs(RowsA);

	vector<vector<byte> > PB;
	MakeMegaProfile_AA(B, PB);

	uint MMLoi, MMLoj, MMLeni, MMLenj;
	string MMPath;
	float MMScore = SWFast_MASM_MegaProf(Mem, MA, PB, g_GapOpen, g_GapExt,
	  MMLoi, MMLoj, MMLeni, MMLenj, MMPath);

	g_PS->m_MASM = &MA;
	g_PS->m_MegaProfile = &PB;

	g_LB = SIZE(B);
	EnumPathsLocal(g_LA, g_LB, OnPath);

	bool AllOk = true;
	++g_N;
	if (g_BestPath != MMPath)
		{
		Log("@PATHDIFF\n");
		++g_NPathDiff;
		AllOk = false;
		}
	if (!feq(g_BestScore, MMScore))
		{
		Log("@SCOREDIFF\n");
		++g_NScoreDiff;
		AllOk = false;
		}
	if (AllOk)
		++g_NAllOk;

	if (!AllOk)
		{
		Log("\n");
		for (uint i = 0; i < g_SeqCountA; ++i)
			Log("%s  >A%u\n", RowsA[i].c_str(), i);
		Log("%s  >B\n", B.c_str());

		Log("  %7.3g  %s  (%u, %u)  Brute\n",
		  g_BestScore, g_BestPath.c_str(), g_BestPosA, g_BestPosB);

		Log("  %7.3g  %s  (%u, %u)  MASM_Mega \n",
		  MMScore, MMPath.c_str(), MMLoi, MMLoj);
		}
	}

static void Test(const string &sA, const string &B)
	{
	ClearBrute();

	vector<string> RowsA;
	Split(sA, RowsA, '|');

	Test(RowsA, B, true);
	}

static void LogPath(const string &sA, const string &B,
  uint PosA, uint PosB, const string &Path)
	{
	Log("___________________________________\n");
	ClearBrute();

	vector<string> RowsA;
	Split(sA, RowsA, '|');

	g_SeqCountA = SIZE(RowsA);
	for (uint i = 0; i < g_SeqCountA; ++i)
		{
		if (i == 0)
			g_LA = SIZE(RowsA[i]);
		else
			asserta(SIZE(RowsA[i]) == g_LA);
		}

	XDPMem Mem;
	MASM &MA = *MakeMASM_AAs(RowsA);
	MA.LogMe();

	vector<vector<byte> > PB;
	MakeMegaProfile_AA(B, PB);
	g_LB = SIZE(B);

	Log("\n");
	for (uint i = 0; i < g_SeqCountA; ++i)
		Log("%s  >A%u\n", RowsA[i].c_str(), i);
	Log("%s  >B\n", B.c_str());
	Log("Path %s\n", Path.c_str());

	g_PS->m_MASM = &MA;
	g_PS->m_MegaProfile = &PB;
	g_PS->m_Trace = true;
	g_PS->m_LA = g_LA;
	g_PS->m_LB = g_LB;
	g_PS->GetLocalScore(PosA, PosB, Path);
	g_PS->m_Trace = false;
	}

static void GetRandSeq(string &s)
	{
	s.clear();
	uint L = g_MinRandL + randu32()%(g_MaxRandL - g_MinRandL + 1);
	asserta(L >= g_MinRandL && L <= g_MaxRandL);
	for (uint i = 0; i < L; ++i)
		s += g_LetterToCharAmino[randu32()%20];
	}

static void GetRandRows(vector<string> &Rows)
	{
	Rows.clear();
	uint L = g_MinRandL + randu32()%(g_MaxRandL - g_MinRandL + 1);
	uint n = g_MinRandN + randu32()%(g_MaxRandN - g_MinRandN + 1);
	asserta(L >= g_MinRandL && L <= g_MaxRandL);
	asserta(n >= g_MinRandN && n <= g_MaxRandN);

	for (uint i = 0; i < n; ++i)
		{
		string Row;
		for (uint j = 0; j < L; ++j)
			{
			if (randu32()%100 < g_RandGapPct)
				Row += '-';
			else
				Row += g_LetterToCharAmino[randu32()%20];
			}
		Rows.push_back(Row);
		}
	for (uint ColIdx = 0; ColIdx < L; ++ColIdx)
		{
		bool AllGaps = true;
		for (uint SeqIdx = 0; SeqIdx < n; ++SeqIdx)
			{
			if (!isgap(Rows[SeqIdx][ColIdx]))
				{
				AllGaps = false;
				break;
				}
			}
		if (AllGaps)
			Rows[randu32()%n][ColIdx] = g_LetterToCharAmino[randu32()%20];
		}
	}

static void TestRandom()
	{
	ClearBrute();
	vector<string> RowsA;
	string B;
	GetRandRows(RowsA);
	GetRandSeq(B);
	Test(RowsA, B, false);
	}

void cmd_test_sw_mm()
	{
	PathScorer_MASM_Mega PS;
	g_PS = &PS;

//WAQHEAW  >A0
//CTFWH  >B
//     5.51  MDDM  (0, 3)  Brute
//     5.25  M  (0, 3)  MASM_Mega 
	Test("WAQHEAW", "CTFWH");
	LogPath("WAQHEAW", "CTFWH", 0, 3, "MDDM");
	return;

	//Test("SEQ|SEQ", "SEQ");
	//Test("SE-|-EQ", "SEQ");
	//Test("SEV-|V-EQ|DDD-", "WEQ");
	//LogPath("SEV-|V-EQ|DDD-", "WEQ", 1, 1, "MDM");
	//LogPath("SEV-|V-EQ|DDD-", "WEQ", 2, 1, "MM");
	//const uint ITERS = 10000;
	//for (uint Iter = 0; Iter < ITERS; ++Iter)
	//	{
	//	ProgressStep(Iter, ITERS, "Running %u ok, %u pathdiff, %u scorediff",
	//	  g_NAllOk, g_NPathDiff, g_NScoreDiff);
	//	TestRandom();
	//	}

	Log("N %u, allok %u, pathdiff %u, scorediff %u\n",
	  g_N, g_NAllOk, g_NPathDiff, g_NScoreDiff);
	}
