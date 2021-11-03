#pragma once

#include "qscorer.h"

class QScorer3
	{
public:
	MSA m_Test1;
	MSA m_Test2;
	MSA m_Ref;
	QScorer m_QS1;
	QScorer m_QS2;

	vector<uint> m_Indexes2;
	const vector<uint> *m_RefCols = 0;
	uint m_RefAlignedColCount = 0;

	vector<pair<uint, uint> > m_Pairs;
	vector<float> m_PairIndexToQ1;
	vector<float> m_PairIndexToQ2;
	vector<float> m_PairIndexToPWC;

public:
	void Run(const string &TestFileName1,
	  const string &TestFileName2,
	  const string &RefFileName,
	  double MaxGapFract);
	void ToTSV(const string &FileName) const;
	void ToTSV(FILE *f) const;

private:
	void TransQ();
	void TransQPair(uint Indexi, uint Indexj);
	};
