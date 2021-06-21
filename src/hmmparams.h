#pragma once

enum HMMTRANS
	{
#define T(x)	HMMTRANS_##x,
#include "hmmtrans.h"
	HMMTRANS_N
	};

#define HMM_ALPHA	"ACDEFGHIKLMNPQRSTVWY"
const uint HMM_ALPHASIZE = 20;

class HMMParams
	{
public:
	bool m_Logs = false;
	vector<float> m_Trans;
	vector<vector<float> > m_Emits;

	vector<string> m_Lines;
	uint m_LineNr = 0;

public:
	HMMParams()
		{
		m_Logs = false;
		m_Trans.resize(HMMTRANS_N, FLT_MAX);
		m_Emits.resize(HMM_ALPHASIZE);
		for (uint i = 0; i < HMM_ALPHASIZE; ++i)
			m_Emits[i].resize(HMM_ALPHASIZE, FLT_MAX);
		}

public:
	void FromProbs(const HMMParams &Probs);
	void FromScores(const HMMParams &Scores);
	void FromStrings(const vector<string> &Lines);
	void FromFile(const string &FileName);
	void FromDefaults();

	void Delta();
	void DeltaEmitProbs(const vector<vector<float> > &Deltas);
	void DeltaTransProbs(float dStartIS, float dStartIL,
	  float dShortOpen, float dShortExtend,
	  float dLongOpen, float dLongExtend);
	void DeltaStartProbs(float dIS, float dIL);
	void DeltaShortGapProbs(float dOpen, float dExtend);
	void DeltaLongGapProbs(float dOpen, float dExtend);

	void AssertProbsValid() const;

	void ToFile(const string &FileName) const;
	void ToPairHMM() const;

	void GetProbs(HMMParams &Probs) const;
	void GetScores(HMMParams &Scores) const;

private:
	bool GetNextLine(string &Line);
	float GetNextProb(const string &Name);

public:
	static void ScoresToProbs(const HMMParams &Scores, HMMParams &Probs);
	static void ProbsToScores(const HMMParams &Probs, HMMParams &Scores);
	static void GetDefaultHMMParams(vector<string> &Lines);
	};
