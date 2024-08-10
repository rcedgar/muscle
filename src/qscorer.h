#pragma once

class QScorer
	{
public:
	const MSA *m_Test;
	const MSA *m_Ref;
	double m_MaxGapFract = 1.0;
	string m_Name;

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
	map<string, uint> m_RefSeqToSeqIndex;
	vector<uint> m_TestColToCount;

	vector<bool> m_RefColIsAligned;
	vector<bool> m_TestColIsAligned;

	double m_RefMSAs_Q = 0;
	uint m_RefMSAs_ComparedColCount = 0;
	vector<uint> m_RefMSAs_TestCols;
	vector<uint> m_RefMSAs_RefCols;
	vector<double> m_RefMSAs_ColQs;

public:
	void Clear();
	bool Run(const string &Name, const MSA &Test, const MSA &Ref);
	void Run(const string &Name, const MultiSequence &Test, const MultiSequence &Ref);
	void CmpRefMSAs(const string &Name, const MSA &Test, const MSA &Ref);

	void InitRefLabels();
	void InitRefToTest();
	void InitRefLabels_bysequence();
	void InitRefToTest_bysequence();
	void InitColPosVecs();
	void InitColPosVecs1(uint i);
	void InitRefCols();
	void InitRefUngappedCounts();
	void DoRefCols();
	void DoRefCol(uint k);
	void SetTestColToBestRefCol();
	void SetTestColIsAligned();
	void SetRefColIsAligned();
	void UpdateRefLetterCounts(vector<vector<uint> > &LetterCountsVec) const;
	void UpdateRefLetterCountsCol(uint k, vector<vector<uint> > &LetterCountsVec) const;

	uint GetRefSeqCount() const { return m_Ref->GetSeqCount(); }
	uint GetTestSeqCount() const { return m_Test->GetSeqCount(); }
	uint GetRefColCount() const { return m_Ref->GetColCount(); }
	uint GetTestColCount() const { return m_Test->GetColCount(); }
	};
