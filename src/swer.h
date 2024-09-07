#pragma once

#include "pathscorer.h"

class SWer
	{
public:
	uint m_LA = UINT_MAX;
	uint m_LB = UINT_MAX;
	string m_A;
	vector<string> m_RowsA;
	string m_B;

public:
	float Run(const string &A, const string &B,
	  uint &LoA, uint &LoB, string &Path);

public:
	virtual float SW(uint &LoA, uint &LoB, string &Path) = 0;
	virtual const char *GetName() const = 0;
	};

class SWer_Enum_Seqs_AA_BLOSUM62 : public SWer
	{
private:
	float m_GapOpen = FLT_MAX;
	float m_GapExt = FLT_MAX;

public:
	float m_BestScore = 0;
	string m_BestPath;
	uint m_BestPosA = UINT_MAX;
	uint m_BestPosB = UINT_MAX;
	PathScorer_AA_BLOSUM62 m_PS;

public:
	virtual float SW(uint &LoA, uint &LoB, string &Path);
	virtual const char *GetName() const { return "Enum_Seqs_AA_BLOSUM62"; };

public:
	void SetGaps(float Open, float Ext);
	};

class SWer_Fast_Seqs_AA_BLOSUM62: public SWer
	{
public:
	float m_GapOpen = FLT_MAX;
	float m_GapExt = FLT_MAX;

public:
	virtual float SW(uint &LoA, uint &LoB, string &Path);
	virtual const char *GetName() const { return "Fast_Seqs_AA_BLOSUM62"; };
	};

class SWer_Simple_Seqs_AA_BLOSUM62: public SWer
	{
public:
	float m_GapOpen = FLT_MAX;
	float m_GapExt = FLT_MAX;
	PathScorer_AA_BLOSUM62 m_PS;

public:
	virtual float SW(uint &LoA, uint &LoB, string &Path);
	virtual const char *GetName() const { return "Simple_Seqs_AA_BLOSUM62"; };
	};

class SWer_Mega_Prof_Seqs: public SWer
	{
public:
	float m_GapOpen = FLT_MAX;
	float m_GapExt = FLT_MAX;

public:
	virtual float SW(uint &LoA, uint &LoB, string &Path);
	virtual const char *GetName() const { return "SWer_Mega_Prof_Seqs"; };
	};
