#pragma once

class QScorer
	{
public:
	const MSA *m_Test;
	const MSA *m_Ref;
	double m_MaxGapFract = 1.0;

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

	float m_Q = 0;
	float m_TC = 0;

	vector<string> m_RefLabels;
	map<string, uint> m_RefLabelToSeqIndex;
	vector<uint> m_TestColToCount;

public:
	void Clear();
	void Run(const MSA &Test, const MSA &Ref);

	void InitRefLabels();
	void InitRefToTest();
	void InitColPosVecs();
	void InitColPosVecs1(uint i);
	void InitRefCols();
	void InitRefUngappedCounts();
	void DoRefCols();
	void DoRefCol(uint k);
	void SetTestColToBestRefCol();
	void UpdateRefLetterCounts(vector<vector<uint> > &LetterCountsVec) const;
	void UpdateRefLetterCountsCol(uint k, vector<vector<uint> > &LetterCountsVec) const;

	uint GetRefSeqCount() const { return m_Ref->GetSeqCount(); }
	uint GetTestSeqCount() const { return m_Test->GetSeqCount(); }
	uint GetRefColCount() const { return m_Ref->GetColCount(); }
	uint GetTestColCount() const { return m_Test->GetColCount(); }
	};
