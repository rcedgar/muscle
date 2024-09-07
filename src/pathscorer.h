#pragma once

#include "masm.h"
#include "mega.h"

class PathScorer
	{
public:
	bool m_Trace = false;
	uint m_LA = UINT_MAX;
	uint m_LB = UINT_MAX;

public:
	float GetLocalScore(uint PosA, uint PosB, const string &Path);
	float GetScore(char FromState, char ToState,
	  uint PosA, uint PosB);
	void Trace(bool On) { m_Trace = On; }

	uint GetLA() const { asserta(m_LA != UINT_MAX); return m_LA; }
	uint GetLB() const { asserta(m_LB != UINT_MAX); return m_LB; }

public:
	virtual float GetMatchScore(uint PosA, uint PosB) = 0;

	virtual float GetScoreMM(uint PosA, uint PosB) = 0;
	virtual float GetScoreMD(uint PosA, uint PosB) = 0;
	virtual float GetScoreMI(uint PosA, uint PosB) = 0;

	virtual float GetScoreDM(uint PosA, uint PosB) = 0;
	virtual float GetScoreDD(uint PosA, uint PosB) = 0;

	virtual float GetScoreIM(uint PosA, uint PosB) = 0;
	virtual float GetScoreII(uint PosA, uint PosB) = 0;
	};

class PathScorer_MASM_Mega : public PathScorer
	{
public:
	const MASM *m_MASM = 0;
	const vector<vector<byte> > *m_MegaProfile = 0;

public:
	virtual float GetMatchScore(uint PosA, uint PosB);
	void Init(MASM &MA, const vector<vector<byte> > &PB);

public:
	virtual float GetScoreMM(uint PosA, uint PosB);
	virtual float GetScoreMD(uint PosA, uint PosB);
	virtual float GetScoreMI(uint PosA, uint PosB);

	virtual float GetScoreDM(uint PosA, uint PosB);
	virtual float GetScoreDD(uint PosA, uint PosB);

	virtual float GetScoreIM(uint PosA, uint PosB);
	virtual float GetScoreII(uint PosA, uint PosB);
	};

class PathScorer_AA_BLOSUM62 : public PathScorer
	{
public:
	float m_GapOpen = FLT_MAX;
	float m_GapExt = FLT_MAX;
	string m_SeqA;
	string m_SeqB;

public:
	float GetMatchScore(uint PosA, uint PosB);

public:
	virtual float GetScoreMM(uint PosA, uint PosB);
	virtual float GetScoreMD(uint PosA, uint PosB);
	virtual float GetScoreMI(uint PosA, uint PosB);

	virtual float GetScoreDM(uint PosA, uint PosB);
	virtual float GetScoreDD(uint PosA, uint PosB);

	virtual float GetScoreIM(uint PosA, uint PosB);
	virtual float GetScoreII(uint PosA, uint PosB);
	};

float SWSimple(PathScorer &PS, uint &LoA, uint &LoB, string &Path);
