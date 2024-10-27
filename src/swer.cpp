#include "myutils.h"
#include "swer.h"
#include "enumpaths.h"
#include "pathscorer.h"
#include "xdpmem.h"
#include "alpha.h"

static void LogM(const char *Msg, const vector<vector<float> > &M)
	{
	const uint LA = SIZE(M);
	asserta(LA > 0);
	const uint LB = SIZE(M[0]);
	Log("LogM(%s)\n", Msg);
	Log("      ");
	for (uint j = 0; j < LB; ++j)
		Log(" %10u", j);
	Log("\n");
	for (uint i = 0; i < LA; ++i)
		{
		Log("%3u | ", i);
		for (uint j = 0; j < LB; ++j)
			{
			float m = M[i][j];
			if (m == FLT_MAX)
				Log(" %10.10s", "*");
			else
				Log(" %10.3g", m);
			}
		Log("\n");
		}
	}

static void CmpFwdM(
  const vector<vector<float> > &M1,
  const vector<vector<float> > &M2)
	{
	const uint LA = SIZE(M1);
	asserta(LA > 0);
	asserta(SIZE(M2) == LA);
	const uint LB = SIZE(M1[0]);
	asserta(SIZE(M2[0]) == LB);
	for (uint i = 0; i < LA; ++i)
		{
		for (uint j = 0; j < LB; ++j)
			{
			float m1 = M1[i][j];
			float m2 = M2[i][j];
			if (!feq(m1, m2))
				{
				Log("i=%u j=%u M1=%.3g M2=%.3g\n", i, j, m1, m2);
				LogM("M1", M1);
				LogM("M2", M2);
				Die("CmpFwdM");
				}
			}
		}
	}

static MASM *MakeMASM_Seq(const string &Seq,
  float GapOpen, float GapExt)
	{
	MultiSequence *Aln = new MultiSequence;
	Sequence *s = NewSequence();
	s->FromString("LABEL", Seq);
	Aln->AddSequence(s, true);
	MASM *M = new MASM;
	Mega::FromMSA_AAOnly(*Aln, GapOpen, GapExt);
	M->FromMSA(*Aln, "MSA", -GapOpen, -GapExt);
	return M;
	}

static MASM *MakeMASM_Rows(const vector<string> &Rows,
  float GapOpen, float GapExt)
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
	Mega::FromMSA_AAOnly(*Aln, GapOpen, GapExt);
	M->FromMSA(*Aln, "Rows", -GapOpen, -GapExt);
	return M;
	}

static void MakeMegaProfile(const string &Seq,
  vector<vector<byte> > &Prof)
	{
	const uint L = SIZE(Seq);
	Prof.clear();
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		vector<byte> Col;
		char c = Seq[Pos];
		byte Letter = g_CharToLetterAmino[c];
		if (Letter >= 20)
			Letter = 0;
		Col.push_back(Letter);
		Prof.push_back(Col);
		}
	}

uint SWer::GetNA(const string &Path) const
	{
	uint n = 0;
	for (uint i = 0; i < SIZE(Path); ++i)
		{
		char c = Path[i];
		if (c == 'M' || c == 'D')
			++n;
		}
	return n;
	}

uint SWer::GetNB(const string &Path) const
	{
	uint n = 0;
	for (uint i = 0; i < SIZE(Path); ++i)
		{
		char c = Path[i];
		if (c == 'M' || c == 'I')
			++n;
		}
	return n;
	}

float SWer::Run(const string &A, const string &B,
  uint &LoA, uint &LoB, string &Path)
	{
	m_RowsA.clear();
	m_A = A;
	m_B = B;
	Split(A, m_RowsA, '|');
	m_LA = SIZE(m_RowsA[0]);
	m_LB = SIZE(m_B);
	float Score = SW(LoA, LoB, Path);
	if (Score <= 0)
		return 0;
	uint NA = GetNA(Path);
	uint NB = GetNB(Path);
	uint HiA = LoA + NA - 1;
	uint HiB = LoB + NB - 1;
	asserta(HiA < m_LA);
	asserta(HiB < m_LB);
	return Score;
	}

static SWer_Enum_Seqs_AA_BLOSUM62 *g_ptrSWer_Brute;
static void OnPath(uint PosA, uint PosB, const string &Path)
	{
	PathScorer &PS = g_ptrSWer_Brute->m_PS;
	PS.m_LA = g_ptrSWer_Brute->m_LA;
	PS.m_LB  = g_ptrSWer_Brute->m_LB;
	float Score = PS.GetLocalScore(PosA, PosB, Path);
	if (Score > g_ptrSWer_Brute->m_BestScore)
		{
		g_ptrSWer_Brute->m_BestScore = Score;
		g_ptrSWer_Brute->m_BestPath = Path;
		g_ptrSWer_Brute->m_BestPosA = PosA;
		g_ptrSWer_Brute->m_BestPosB = PosB;
		}
	}

void SWer_Enum_Seqs_AA_BLOSUM62::SetGaps(float Open, float Ext)
	{
	m_GapOpen = Open;
	m_GapExt = Ext;
	m_PS.m_GapOpen = Open;
	m_PS.m_GapExt = Ext;
	}

float SWer_Enum_Seqs_AA_BLOSUM62::SW(uint &LoA, uint &LoB, string &Path)
	{
	m_BestScore = 0;
	m_BestPath.clear();
	m_BestPosA = UINT_MAX;
	m_BestPosB = UINT_MAX;

	g_ptrSWer_Brute = this;
	g_ptrSWer_Brute->m_PS.m_GapOpen = m_GapOpen;
	g_ptrSWer_Brute->m_PS.m_SeqA = m_A;
	g_ptrSWer_Brute->m_PS.m_SeqB = m_B;
	EnumPathsLocal(m_LA, m_LB, OnPath);
	LoA = m_BestPosA;
	LoB = m_BestPosB;
	Path = m_BestPath;
	return m_BestScore;
	}

float SWer_Fast_Seqs_AA_BLOSUM62::SW(uint &LoA, uint &LoB, string &Path)
	{
	float SWFast_Strings_BLOSUM62(XDPMem &Mem,
	  const string &A, const string &B, float Open, float Ext,
	  uint &Loi, uint &Loj, uint &Leni, uint &Lenj, string &Path);
	asserta(m_GapOpen != FLT_MAX && m_GapOpen < 0);
	asserta(m_GapExt != FLT_MAX && m_GapExt < 0);

	m_PS.m_GapOpen = m_GapOpen;
	m_PS.m_GapExt = m_GapExt;
	m_PS.m_LA = m_LA;
	m_PS.m_LB = m_LB;
	m_PS.m_SeqA = m_A;
	m_PS.m_SeqB = m_B;

	XDPMem Mem;
	uint Leni, Lenj;
	float Score = SWFast_Strings_BLOSUM62(Mem, m_A, m_B, m_GapOpen, m_GapExt,
	  LoA, LoB, Leni, Lenj, Path);
	return Score;
	}

float SWer_Simple_Seqs_AA_BLOSUM62::SW(uint &LoA, uint &LoB, string &Path)
	{
	float SWSimple(PathScorer &PS, uint &LoA, uint &LoB, string &Path);
	asserta(m_GapOpen != FLT_MAX && m_GapOpen < 0);
	asserta(m_GapExt != FLT_MAX && m_GapExt < 0);

	m_PS.m_GapOpen = m_GapOpen;
	m_PS.m_GapExt = m_GapExt;
	m_PS.m_SeqA = m_A;
	m_PS.m_SeqB = m_B;
	m_PS.m_LA = SIZE(m_A);
	m_PS.m_LB = SIZE(m_B);

	float Score = SWSimple(m_PS, LoA, LoB, Path);
	return Score;
	}

float SWer_MASM_Mega_Seqs::SW(uint &LoA, uint &LoB, string &Path)
	{
	float SWFast_MASM_MegaProf(XDPMem &Mem, const MASM &MA,
	  const vector<vector<byte> > &PB, float Open, float Ext,
	  uint &Loi, uint &Loj, uint &Leni, uint &Lenj, string &Path);

	asserta(m_GapOpen != FLT_MAX && m_GapOpen < 0);
	asserta(m_GapExt != FLT_MAX && m_GapExt < 0);

	MASM *MA = MakeMASM_Seq(m_A, m_GapOpen, m_GapExt);
	vector<vector<byte> > &PB = *new vector<vector<byte> >;
	MakeMegaProfile(m_B, PB);

	m_PS.m_LA = m_LA;
	m_PS.m_LB = m_LB;
	m_PS.m_MASM = MA;
	m_PS.m_MegaProfile = &PB;

	XDPMem Mem;
	uint Leni, Lenj;
	float Score = SWFast_MASM_MegaProf(Mem, *MA, PB, m_GapOpen, m_GapExt,
	  LoA, LoB, Leni, Lenj, Path);

	return Score;
	}

float SWer_Simple_MASM_Mega::SW(uint &LoA, uint &LoB, string &Path)
	{
	float SWSimple(PathScorer &PS, uint &LoA, uint &LoB, string &Path);
	float SWEnumDPFwdM(PathScorer &PS, uint LA, uint LB, uint &LoA, uint &LoB,
	  string &Path, vector<vector<float> > &FwdM);
	float SWSimpleFwdM(PathScorer &PS, uint &LoA, uint &LoB, string &Path,
	  vector<vector<float> > &FwdM);

	asserta(m_GapOpen != FLT_MAX && m_GapOpen < 0);
	asserta(m_GapExt != FLT_MAX && m_GapExt < 0);

	vector<vector<byte> > &PB = *new vector<vector<byte> >;
	MakeMegaProfile(m_B, PB);
	m_PS.m_MASM = MakeMASM_Rows(m_RowsA, m_GapOpen, m_GapExt);
	m_PS.m_MegaProfile = &PB;
	m_PS.m_LA = SIZE(m_RowsA[0]);
	m_PS.m_LB = SIZE(m_B);

	float Score = SWSimple(m_PS, LoA, LoB, Path);

	vector<vector<float> > EnumFwdM;
	vector<vector<float> > SimpleFwdM;
	string Path2;
	string Path3;
	float Score2 = SWEnumDPFwdM(m_PS, m_LA, m_LB, LoA, LoB, Path2, EnumFwdM);
	float Score3 = SWSimpleFwdM(m_PS, LoA, LoB, Path3, SimpleFwdM);
	//Log("Score %.3g, Score2 %.3g, Score3 %.3g\n", Score, Score2, Score3);
	asserta(feq(Score, Score2));
	asserta(feq(Score, Score3));
	//CmpFwdM(EnumFwdM, SimpleFwdM);

	return Score;
	}

float SWer_MASM_Mega::SW(uint &LoA, uint &LoB, string &Path)
	{
	float SWFast_MASM_MegaProf(XDPMem &Mem, const MASM &MA,
	  const vector<vector<byte> > &PB, float Open, float Ext,
	  uint &Loi, uint &Loj, uint &Leni, uint &Lenj, string &Path);

	asserta(m_GapOpen != FLT_MAX && m_GapOpen < 0);
	asserta(m_GapExt != FLT_MAX && m_GapExt < 0);

	MASM *MA = MakeMASM_Rows(m_RowsA, m_GapOpen, m_GapExt);
	vector<vector<byte> > PB;
	MakeMegaProfile(m_B, PB);

	XDPMem Mem;
	uint Leni, Lenj;
	float Score = SWFast_MASM_MegaProf(Mem, *MA, PB, m_GapOpen, m_GapExt,
	  LoA, LoB, Leni, Lenj, Path);

	return Score;
	}

//static SWer_Enum_MASM_Mega *g_ptrSWer_Brute2;
//static void OnPath2(uint PosA, uint PosB, const string &Path)
//	{
//	PathScorer &PS = g_ptrSWer_Brute2->m_PS;
//	PS.m_LA = g_ptrSWer_Brute2->m_LA;
//	PS.m_LB  = g_ptrSWer_Brute2->m_LB;
//	float Score = PS.GetLocalScore(PosA, PosB, Path);
//	if (Score > g_ptrSWer_Brute2->m_BestScore)
//		{
//		g_ptrSWer_Brute2->m_BestScore = Score;
//		g_ptrSWer_Brute2->m_BestPath = Path;
//		g_ptrSWer_Brute2->m_BestPosA = PosA;
//		g_ptrSWer_Brute2->m_BestPosB = PosB;
//		}
//	}

float SWer_Enum_MASM_Mega::SW(uint &LoA, uint &LoB, string &Path)
	{
	asserta(m_GapOpen != FLT_MAX && m_GapOpen < 0);
	asserta(m_GapExt != FLT_MAX && m_GapExt < 0);

	vector<string> RowsA;
	Split(m_A, RowsA, '|');
	MASM *MA = MakeMASM_Rows(RowsA, m_GapOpen, m_GapExt);
	vector<vector<byte> > &PB = *new vector<vector<byte> >;
	MakeMegaProfile(m_B, PB);

	m_PS.m_LA = m_LA;
	m_PS.m_LB = m_LB;
	m_PS.m_MASM = MA;
	m_PS.m_MegaProfile = &PB;

	float SWEnumDP(PathScorer &PS, uint LA, uint LB, uint &LoA, uint &LoB,
	  string &Path);

	float Score = SWEnumDP(m_PS, m_LA, m_LB, LoA, LoB, Path);
	return Score;
	}

float SWer_PS::SW(uint &LoA, uint &LoB, string &Path)
	{
	float SWPS(XDPMem &Mem, PathScorer &PS, uint &Loi, uint &Loj, string &Path);
	float SWSimple2(XDPMem &Mem, PathScorer &PS, uint &LoA, uint &LoB, string &Path);

	asserta(m_PS != 0);
	PathScorer_AA_BLOSUM62 *PS = (PathScorer_AA_BLOSUM62 *) m_PS;
	PS->m_SeqA = m_A;
	PS->m_SeqB = m_B;
	PS->m_LA = m_LA;
	PS->m_LB = m_LB;

	XDPMem Mem;
	float Score = SWPS(Mem, *m_PS, LoA, LoB, Path);

	XDPMem Mem2;
	uint LoA2, LoB2;
	string Path2;
	float Score2 = SWSimple2(Mem2, *m_PS, LoA2, LoB2, Path);
	asserta(feq(Score, Score2));
	return Score;
	}
