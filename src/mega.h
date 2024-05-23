#pragma once

class Mega
	{
public:
	string m_FileName;
	vector<string> m_Lines;
	vector<string> m_FeatureNames;
	vector<float> m_Weights;
	vector<uint> m_AlphaSizes;
	vector<vector<vector<float> > > m_SubstMxVec;
	vector<vector<float> > m_FreqsVec;
	vector<string> m_Labels;
	vector<vector<vector<byte> > > m_Profiles;
	vector<string> m_Seqs;
	uint m_NextLineNr = 0;

public:
	void FromFile(const string &FileName);
	const string &GetNextLine();
	void GetNextFields(vector<string> &Fields,
	  uint ExpectedNrFields = UINT_MAX);
	};
