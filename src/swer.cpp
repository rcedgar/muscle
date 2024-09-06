#include "myutils.h"
#include "swer.h"
#include "enumpaths.h"
#include "pathscorer.h"
#include "xdpmem.h"

static SWer_Enum_Seqs_AA_BLOSUM62 *g_ptrSWer_Brute;

float SWer::Run(const string &A, const string &B,
  uint &LoA, uint &LoB, string &Path)
	{
	m_RowsA.clear();
	m_A = A;
	m_B = B;
	m_LA = SIZE(m_A);
	m_LB = SIZE(m_B);
	size_t n = A.find('|');
	if (n == string::npos)
		Split(A, m_RowsA, '|');
	float Score = SW(LoA, LoB, Path);
	return Score;
	}

static void OnPath(uint PosA, uint PosB, const string &Path)
	{
	PathScorer &PS = g_ptrSWer_Brute->m_PS;
	uint LA = g_ptrSWer_Brute->m_LA;
	uint LB = g_ptrSWer_Brute->m_LB;
	float Score = PS.GetLocalScore(PosA, PosB, LA, LB, Path);
	if (Score > g_ptrSWer_Brute->m_BestScore)
		{
		g_ptrSWer_Brute->m_BestScore = Score;
		g_ptrSWer_Brute->m_BestPath = Path;
		g_ptrSWer_Brute->m_BestPosA = PosA;
		g_ptrSWer_Brute->m_BestPosB = PosB;
		}
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
	return m_BestScore;
	}

float SWer_Fast_Seqs_AA_BLOSUM62::SW(uint &LoA, uint &LoB, string &Path)
	{
	float SWFast_Strings_BLOSUM62(XDPMem &Mem,
	  const string &A, const string &B, float Open, float Ext,
	  uint &Loi, uint &Loj, uint &Leni, uint &Lenj, string &Path);
	asserta(m_GapOpen != FLT_MAX && m_GapOpen < 0);
	asserta(m_GapExt != FLT_MAX && m_GapExt < 0);

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

	float Score = SWSimple(m_PS, LoA, LoB, Path);
	return Score;
	}
