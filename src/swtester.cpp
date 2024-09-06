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
	if (X_Score != Y_Score)
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
	Log("%s\n", Msg);
	Log("A: %s\n", m_A.c_str());
	Log("B: %s\n", m_B.c_str());
	Log("  %.3g/%.3g", X_Score, Y_Score);
	Log("  %u,%u", X_LoA, Y_LoA);
	Log("  %s,%s", X_Path.c_str(), Y_Path.c_str());
	Log("\n");
	}

void SWTester::GetRandomSeq(uint L, string &s)
	{
	s.clear();
	for (uint i = 0; i < L; ++i)
		s += g_LetterToCharAmino[randu32()%20];
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
