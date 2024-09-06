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
	void SetX(SWer &X) { m_X = &X; }
	void SetY(SWer &Y) { m_Y = &Y; }
	void RunX(const string &A, const string &B);
	void RunY(const string &A, const string &B);
	void RunXY(SWer &X, SWer &Y, const string &A, const string &B);
	void RunAB(const string &A, const string &B);
	void CmpXY();
	void LogResult(const char *Msg) const;
	void RunRandomSeqs(uint MinL, uint MaxL);
	void RunRandomMSASeq(uint MinN, int MaxN, uint MinL, uint MaxL);
	void GetRandomSeq(uint L, string &Seq);
	};
