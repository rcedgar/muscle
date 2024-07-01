#pragma once

class QScorer
	{
public:
	string m_Name;
	const MSA *m_Test = 0;
	const MSA *m_Ref = 0;
	double m_MaxGapFract = 1.0;
	bool m_MissingTestSeqOk = false;

	uint m_RefAlignedColCount = 0;

	vector<string> m_Labels;
	vector<uint> m_RefSeqIndexes;
	vector<uint> m_TestSeqIndexes;
	vector<uint> m_RefSeqIndexToTestSeqIndex;
	vector<uint> m_TestSeqIndexToRefSeqIndex;
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

        map<string, uint> m_RefSeqToSeqIndex;

public:
	QScorer()
		{
		m_MissingTestSeqOk = opt(missingtestseqok);
		}

	void Clear();
	void Run(const string &Name, const MSA &Test, const MSA &Ref);
	void Run(const string &Name, const MultiSequence &Test, const MultiSequence &Ref);

	void InitRefLabels();
	void InitRefToTest();
	void InitRefToTest_BySequence();
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
private:
	map<string, uint>::const_iterator _FullScanRefSeqs(const string & TestSeq);
	};
