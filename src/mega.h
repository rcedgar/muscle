#pragma once

#include <unordered_map>

class Mega
	{
public:
	static string m_FileName;
	static vector<string> m_Lines;
	static vector<string> m_FeatureNames;
	static vector<float> m_Weights;
	static vector<uint> m_AlphaSizes;
	static unordered_map<string, uint> m_LabelToIdx;

// log(P_i) for each letter (for HMM Insert states)
	static vector<vector<float> > m_LogProbsVec;

// log(P_ij) for each letter pair (for HMM Match state)
	static vector<vector<vector<float> > > m_LogProbMxVec;

	static vector<string> m_Labels;
	static vector<vector<vector<byte> > > m_Profiles;
	static vector<string> m_Seqs;
	static uint m_NextLineNr;
	static uint m_FeatureCount;

public:
	void FromFile(const string &FileName);
	uint GetFeatureCount() const { return m_FeatureCount; }
	const string &GetNextLine();
	void GetNextFields(vector<string> &Fields,
	  uint ExpectedNrFields = UINT_MAX);
	float GetInsScore(const vector<vector<byte> > &Profile, uint Pos) const;
	float GetMatchScore(
	  const vector<vector<byte> > &ProfileX, uint PosX,
	  const vector<vector<byte> > &ProfileY, uint PosY) const;
	void CalcLogProbsMx(const vector<vector<float > > &FreqsMx,
	  vector<vector<float > > &LogProbMx) const;
	void CalcMarginalFreqs(const vector<vector<float > > &FreqsMx,
	  vector<float> &Freqs) const;
	void LogFeatureParams(uint Idx) const;
	void LogMx(const string &Name, const vector<vector<float> > &Mx) const;
	void LogVec(const string &Name, const vector<float> &Vec) const;
	void AssertSymmetrical(const vector<vector<float> > &Mx) const;

public:
	static uint GetGSIByLabel(const string &Label);
	static const string &GetLabelByGSI(uint GSI);
	static const vector<vector<byte> > *GetProfileByGSI(uint GSI);
	static const vector<vector<byte> > *GetProfileByLabel(const string &Label);
	};
