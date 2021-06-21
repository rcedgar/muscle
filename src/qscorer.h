#pragma once

class QScorer
	{
public:
	string m_TestName;
	string m_RefName;

	const MSA *m_Test;
	const MSA *m_Ref;
	uint m_RefAlignedColCount = 0;

	vector<string> m_Labels;
	vector<uint> m_RefSeqIndexes;
	vector<uint> m_TestSeqIndexes;
	vector<uint> m_RefSeqIndexToTestSeqIndex;
	vector<uint> m_RefCols;
	vector<uint> m_RefUngappedCounts;

	vector<vector<uint> > m_PosToTestColVec;
	vector<vector<uint> > m_PosToRefColVec;

	vector<vector<uint> > m_TestColToPosVec;
	vector<vector<uint> > m_RefColToPosVec;

	vector<vector<uint> > m_RefColToTestColVec;
	vector<uint> m_TestColToBestRefCol;
	vector<float> m_MaxFracts;
	vector<uint> m_BestTestCols;

	uint64 m_TotalPairs = 0;
	uint64 m_TotalCols = 0;

	uint64 m_CorrectPairs = 0;
	uint m_CorrectCols = 0;
	uint m_CorrectCols90 = 0;

	float m_Q = 0;
	float m_TC = 0;
	float m_TC90 = 0;
	float m_CF = 0;
	float m_Acc = 0;

public:
	void Run(const string &TestName, const string &RefName,
	  const MSA &Test, const MSA &Ref);
	void Report(FILE *f) const;
	void Report1(FILE *f, uint Index) const;
	};
