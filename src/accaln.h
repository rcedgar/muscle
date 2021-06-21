#pragma once

class AccAln
	{
public:
	vector<string> m_Labels1;
	vector<string> m_Labels2;
	vector<string> m_Rows1;
	vector<string> m_Rows2;
	vector<string> m_EADigits1;
	vector<string> m_EADigits2;
	vector<vector<float> > m_EAs;
	FILE *m_f = 0;
	uint m_LineNr = 0;
	uint m_PairIndex = 0;
	string m_FileName;
	string m_RefFileName;
	MSA m_RefMSA;
	map<string, uint> m_RefLabelToIndex;
	vector<string> m_RefRows1;
	vector<string> m_RefRows2;
	vector<vector<bool> > m_TestColToCorrectVec;
	vector<vector<bool> > m_TestColToAlignedVec;

public:
	void Clear()
		{
		m_Labels1.clear();
		m_Labels2.clear();
		m_Rows1.clear();
		m_Rows2.clear();
		m_EAs.clear();
		m_EADigits1.clear();
		m_EADigits2.clear();
		m_LineNr = 0;
		m_FileName.clear();
		m_RefFileName.clear();
		m_RefMSA.Clear();
		m_RefLabelToIndex.clear();
		m_RefRows1.clear();
		m_RefRows2.clear();
		m_TestColToCorrectVec.clear();
		m_TestColToAlignedVec.clear();
		}

	void FromFile(const string &FileName);
	void ReadRef(const string &FileName);
	uint GetPairCount() const;
	bool NextPair();
	bool GetBlankLine();
	void GetAlnRow(string &Row);
	void AlnRef();
	void AlnRefPair(uint PairIndex);
	void GetCorrect(uint PairIndex, vector<bool> &ColToCorrect,
	  vector<bool> &ColToAligned);
	void WriteRefAln(FILE *f, uint PairIndex) const;
	void WriteRefAln(const string &FileName) const;
	float GetAlignedEA(uint PairIndex) const;
	float GetQ(uint PairIndex) const;
	};
