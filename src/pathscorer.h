#pragma once

#include "masm.h"
#include "mega.h"

class PathScorer
	{
public:
	float GetGlobalScore(uint LA, uint LB, const string &Path);
	float GetLocalScore(uint PosA, uint PosB, uint LA, uint LB, 
	  const string &Path);
	float GetScore(char FromState, char ToState,
	  uint PosA, uint PosB);
	void TermizePath(string &Path);

public:
	virtual float GetScoreSM(uint PosA, uint PosB) = 0;
	virtual float GetScoreSd(uint PosA, uint PosB) = 0;
	virtual float GetScoreSi(uint PosA, uint PosB) = 0;

	virtual float GetScoreMM(uint PosA, uint PosB) = 0;
	virtual float GetScoreMD(uint PosA, uint PosB) = 0;
	virtual float GetScoreMd(uint PosA, uint PosB) = 0;
	virtual float GetScoreMI(uint PosA, uint PosB) = 0;
	virtual float GetScoreMi(uint PosA, uint PosB) = 0;

	virtual float GetScoreDM(uint PosA, uint PosB) = 0;
	virtual float GetScoredM(uint PosA, uint PosB) = 0;
	virtual float GetScoreDD(uint PosA, uint PosB) = 0;
	virtual float GetScoredd(uint PosA, uint PosB) = 0;

	virtual float GetScoreIM(uint PosA, uint PosB) = 0;
	virtual float GetScoreiM(uint PosA, uint PosB) = 0;
	virtual float GetScoreII(uint PosA, uint PosB) = 0;
	virtual float GetScoreii(uint PosA, uint PosB) = 0;

	virtual float GetScoreME(uint PosA, uint PosB) = 0;
	virtual float GetScoredE(uint PosA, uint PosB) = 0;
	virtual float GetScoreiE(uint PosA, uint PosB) = 0;
	};

class PathScorer_MASM_Mega
	{
public:
	MASM *m_MASM = 0;
	const vector<vector<byte> > *m_MegaProfile = 0;
	float GetMatchScore(uint PosA, uint PosB);

public:
	virtual float GetScoreSM(uint PosA, uint PosB);
	virtual float GetScoreSd(uint PosA, uint PosB);
	virtual float GetScoreSi(uint PosA, uint PosB);

	virtual float GetScoreMM(uint PosA, uint PosB);
	virtual float GetScoreMD(uint PosA, uint PosB);
	virtual float GetScoreMd(uint PosA, uint PosB);
	virtual float GetScoreMI(uint PosA, uint PosB);
	virtual float GetScoreMi(uint PosA, uint PosB);

	virtual float GetScoreDM(uint PosA, uint PosB);
	virtual float GetScoredM(uint PosA, uint PosB);
	virtual float GetScoreDD(uint PosA, uint PosB);
	virtual float GetScoredd(uint PosA, uint PosB);

	virtual float GetScoreIM(uint PosA, uint PosB);
	virtual float GetScoreiM(uint PosA, uint PosB);
	virtual float GetScoreII(uint PosA, uint PosB);
	virtual float GetScoreii(uint PosA, uint PosB);

	virtual float GetScoreME(uint PosA, uint PosB);
	virtual float GetScoredE(uint PosA, uint PosB);
	virtual float GetScoreiE(uint PosA, uint PosB);
	};
