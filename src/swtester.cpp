#include "myutils.h"
#include "swtester.h"
#include "alpha.h"

void SWTester::RunX(const string &A, const string &B)
	{
	m_A = A;
	m_B = B;
	X_Score = m_X->Run(m_A, m_B, X_LoA, X_LoB, X_Path);
	}

void SWTester::RunY(const string &A, const string &B)
	{
	m_A = A;
	m_B = B;
	Y_Score = m_Y->Run(m_A, m_B, Y_LoA, Y_LoB, Y_Path);
	}

void SWTester::RunXY(SWer &X, SWer &Y, const string &A, const string &B)
	{
	SetX(X);
	SetY(Y);
	RunX(A, B);
	RunY(A, B);
	CmpXY();
	}

void SWTester::RunXAB(SWer &X, const string &A, const string &B, bool Trace)
	{
	SetX(X);
	RunX(A, B);
	PathScorer &PS = *m_X->GetPS();
	PS.m_Trace = Trace;
	float ScorePS = PS.GetLocalScore(X_LoA, X_LoB, X_Path);
	Log("\nRunXAB(%s)\n", m_X->GetName());
	Log("A %s\n", A.c_str());
	Log("B %s\n", B.c_str());
	Log("ScoreSW %.3g, ScorePS %.3g", X_Score, ScorePS);
	Log("   %s\n", X_Path.c_str());
	}

void SWTester::RunAB( const string &A, const string &B)
	{
	RunX(A, B);
	RunY(A, B);
	CmpXY();
	}

void SWTester::CmpXY()
	{
	bool Agree = true;
	++m_N;
	if (X_Score == 0 && Y_Score == 0)
		{
		++m_NAgree;
		return;
		}

	if (!feq(X_Score, Y_Score))
		{
		if (Agree)
			LogResult("@SCOREDIFF");
		Agree = false;
		++m_NScoreDiff;
		}
	if (X_Path != Y_Path)
		{
		if (Agree)
			LogResult("@PATHDIFF");
		Agree = false;
		++m_NPathDiff;
		}
	if (X_LoA != Y_LoA || X_LoB != Y_LoB)
		{
		if (Agree)
			LogResult("@POSDIFF");
		Agree = false;
		++m_NPosDiff;
		}
	if (Agree)
		++m_NAgree;
	}

void SWTester::LogResult(const char *Msg) const
	{
	Log("%s %s %s\n", Msg, m_X->GetName(), m_Y->GetName());
	Log("A: %s\n", m_A.c_str());
	Log("B: %s\n", m_B.c_str());
	Log("  %.3g/%.3g", X_Score, Y_Score);
	Log("  loa %u,%u", X_LoA, Y_LoA);
	Log("  lob %u,%u", X_LoB, Y_LoB);
	Log("  %s,%s", X_Path.c_str(), Y_Path.c_str());
	Log("\n");
	}

void SWTester::GetRandomSeq(uint L, string &s)
	{
	s.clear();
	for (uint i = 0; i < L; ++i)
		s += g_LetterToCharAmino[randu32()%20];
	}

void SWTester::RunRandomSeqsIters(uint MinL, uint MaxL, uint Iters)
	{
	for (uint Iter = 0; Iter < Iters; ++Iter)
		{
		ProgressStep(Iter, Iters, "RunRandomSeqsIters");
		RunRandomSeqs(MinL, MaxL);
		}
	}

void SWTester::RunRandomSeqs(uint MinL, uint MaxL)
	{
	uint LA = MinL + randu32()%(MaxL - MinL + 1);
	uint LB = MinL + randu32()%(MaxL - MinL + 1);
	string A, B;
	GetRandomSeq(LA, A);
	GetRandomSeq(LB, B);
	RunAB(A, B);
	}

void SWTester::RunRandomMSASeq(uint MinN, int MaxN, uint MinL, uint MaxL)
	{
	uint N = MinN + randu32()%(MaxN - MinN + 1);
	uint LA = MinL + randu32()%(MaxL - MinL + 1);
	uint LB = MinL + randu32()%(MaxL - MinL + 1);
	string A;
	for (uint i = 0; i < N; ++i)
		{
		if (i > 0)
			A += '|';
		string s;
		GetRandomSeq(LA, s);
		A += s;
		}

	string B;
	GetRandomSeq(LB, B);
	RunAB(A, B);
	}

void SWTester::Stats()
	{
	ProgressLog("\n");
	if (m_X != 0)
		ProgressLog("%10.10s  %s\n", "X", m_X->GetName());
	if (m_Y != 0)
		ProgressLog("%10.10s  %s\n", "Y", m_Y->GetName());
	ProgressLog("%10u  Tests\n", m_N);
	ProgressLog("%10u  Agree\n", m_NAgree);
	ProgressLog("%10u  Score diff\n", m_NScoreDiff);
	ProgressLog("%10u  Path diff\n", m_NPathDiff);
	ProgressLog("%10u  Pos diff\n", m_NPosDiff);
	}
