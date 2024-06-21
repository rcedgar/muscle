#pragma once

class Mega
	{
public:
	string m_FileName;
	vector<string> m_Lines;
	vector<string> m_FeatureNames;
	vector<float> m_Weights;
	vector<uint> m_AlphaSizes;

// log(P_i) for each letter (for HMM Insert states)
	vector<vector<float> > m_LogProbsVec;

// log(P_ij) for each letter pair (for HMM Match state)
	vector<vector<vector<float> > > m_LogProbMxVec;

	vector<string> m_Labels;
	vector<vector<vector<byte> > > m_Profiles;
	vector<string> m_Seqs;
	uint m_NextLineNr = 0;
	uint m_FeatureCount = 0;

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
	};
