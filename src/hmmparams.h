#pragma once

enum HMMTRANS
	{
#define T(x)	HMMTRANS_##x,
#include "hmmtrans.h"
	HMMTRANS_N
	};

static const string AMINO_ALPHA = "ACDEFGHIKLMNPQRSTVWY";
static const string NT_ALPHA = "ACGT";

//#define HMM_ALPHA	"ACDEFGHIKLMNPQRSTVWY"
//const uint HMM_ALPHASIZE = 20;

const float DEFAULT_PERTURB_VAR = 0.25f;

class HMMParams
	{
public:
	bool m_Logs = false;
	uint m_LineNr = 0;
	float m_Var;

	vector<float> m_Trans;
	vector<vector<float> > m_Emits;
	vector<string> m_Lines;
	string m_Alpha;

public:
	HMMParams()
		{
		m_Logs = false;
		m_LineNr = 0;
		m_Var = DEFAULT_PERTURB_VAR;
		}

public:
	void Clear()
		{
		m_Logs = false;
		m_LineNr = 0;
		m_Var = DEFAULT_PERTURB_VAR;

		m_Trans.clear();
		m_Emits.clear();
		m_Lines.clear();
		m_Alpha.clear();
		}

	uint GetAlphaSize() const
		{
		uint n = SIZE(m_Alpha);
		asserta(n == 4 || n == 20);
		return n;
		}

	void FromParams(const HMMParams &Params, bool AsProbs);
	void FromStrings(const vector<string> &Lines);
	void FromFile(const string &FileName);
	void FromDefaults(bool Nucleo);
	void CmdLineUpdate();

	void PerturbProbs(uint Seed);
	void ToSingleAffineProbs(HMMParams &Params);

	void Normalize();
	void NormalizeStart();
	void NormalizeMatch();
	void NormalizeShortGap();
	void NormalizeLongGap();
	void NormalizeEmit();

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
	static void GetDefaultHMMParams(bool Nucleo, vector<string> &Lines);
	static void GetDefaultHMMParams_Amino(vector<string> &Lines);
	static void GetDefaultHMMParams_Nucleo(vector<string> &Lines);
	static void Compare(const HMMParams &HP1, const HMMParams &HP2,
	  float &MeanTransDelta, float &MeanEmitDelta);
	};
