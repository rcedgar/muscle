#pragma once

#include "pathscorer.h"
#include "swer.h"

class SWTester
	{
public:
	SWer *m_X = 0;
	SWer *m_Y = 0;

	uint m_N = 0;
	uint m_NAgree = 0;
	uint m_NScoreDiff = 0;
	uint m_NPathDiff = 0;
	uint m_NPosDiff = 0;
	uint m_NPSScoreOk = 0;
	uint m_NPSScoreDiff = 0;

	string m_A;
	string m_B;

	float X_Score = FLT_MAX;
	float Y_Score = FLT_MAX;

	uint X_LoA = UINT_MAX;
	uint Y_LoA = UINT_MAX;

	uint X_LoB = UINT_MAX;
	uint Y_LoB = UINT_MAX;

	string X_Path;
	string Y_Path;

public:
	void ClearStats()
		{
		m_N = 0;
		m_NAgree = 0;
		m_NScoreDiff = 0;
		m_NPathDiff = 0;
		m_NPosDiff = 0;
		m_NPSScoreOk = 0;
		m_NPSScoreDiff = 0;
		}
	void SetX(SWer &X) { m_X = &X; }
	void SetY(SWer &Y) { m_Y = &Y; }
	void RunXAB(SWer &X, const string &A, const string &B, bool Trace);
	void RunX(const string &A, const string &B);
	void RunY(const string &A, const string &B);
	void RunXY(SWer &X, SWer &Y, const string &A, const string &B);
	void RunAB(const string &A, const string &B);
	void CmpXY();
	void LogResult(const char *Msg) const;
	void RunRandomSeqs(uint MinL, uint MaxL);
	void RunRandomMSASeq(uint MinN, uint MaxN, uint MinL, uint MaxL);
	void GetRandomSeq(uint L, string &Seq, bool WithGaps);
	void RunRandomSeqsIters(uint MinL, uint MaxL, uint Iters);
	void RunRandomMSASeqIters(uint MinN, uint MaxN, uint MinL, uint MaxL, uint Iters);
	void Stats();
	void FixGaps(string &AlnBar);
	};
