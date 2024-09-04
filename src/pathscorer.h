#pragma once

#include "masm.h"
#include "mega.h"

class PathScorer
	{
public:
	float GetLocalScore(uint PosA, uint PosB, uint LA, uint LB, 
	  const string &Path);
	float GetScore(char FromState, char ToState,
	  uint PosA, uint PosB);

public:
	virtual float GetScoreMM(uint PosA, uint PosB) = 0;
	virtual float GetScoreMD(uint PosA, uint PosB) = 0;
	virtual float GetScoreMI(uint PosA, uint PosB) = 0;

	virtual float GetScoreDM(uint PosA, uint PosB) = 0;
	virtual float GetScoreDD(uint PosA, uint PosB) = 0;

	virtual float GetScoreIM(uint PosA, uint PosB) = 0;
	virtual float GetScoreII(uint PosA, uint PosB) = 0;
	};

class PathScorer_MASM_Mega
	{
public:
	MASM *m_MASM = 0;
	const vector<vector<byte> > *m_MegaProfile = 0;
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
