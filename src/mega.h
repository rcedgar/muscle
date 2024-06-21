#pragma once

class Mega
	{
public:
	string m_FileName;
	vector<string> m_Lines;
	vector<string> m_FeatureNames;
	vector<float> m_Weights;
	vector<uint> m_AlphaSizes;
	vector<vector<vector<float> > > m_FreqsMxVec;
	vector<vector<vector<float> > > m_LogProbMxVec;
	vector<vector<float> > m_LogProbsVec;
	vector<string> m_Labels;
	vector<vector<vector<byte> > > m_Profiles;
	vector<string> m_Seqs;
	uint m_NextLineNr = 0;
	uint m_FeatureCount = 0;

public:
	void FromFile(const string &FileName);
	const string &GetNextLine();
	void GetNextFields(vector<string> &Fields,
	  uint ExpectedNrFields = UINT_MAX);
	float GetInsScore(const vector<vector<byte> > &Profile, uint Pos) const;
	float GetMatchScore(
	  const vector<vector<byte> > &ProfileX, uint PosX,
	  const vector<vector<byte> > &ProfileY, uint PosY) const;
	void CalcLogProbsMx(const vector<vector<float > > &FreqsMx,
	  vector<vector<float > > &LogProbMx);
	void CalcMarginalFreqs(const vector<vector<float > > &FreqsMx,
	  vector<float> &Freqs);
	void LogFeatureParams(uint Idx) const;
	void LogMx(const string &Name,
	  const vector<vector<float> > &Mx) const;
	void AssertSymmetrical(const vector<vector<float> > &Mx) const;
	};
