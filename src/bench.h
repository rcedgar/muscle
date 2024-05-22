#pragma once

#include "multisequence.h"
#include "muscle3.h"
#include "qscorer.h"
#include "qscorer2.h"

class Bench
	{
public:
	const M3AlnParams *m_AP = 0;
	vector<string> m_RefNames;
	vector<MultiSequence *> m_Refs;
	vector<MultiSequence *> m_Inputs;
	vector<Muscle3 *> m_M3s;
	vector<QScorer *> m_QSs;
	vector<QScorer2 *> m_QS2s;
	uint m_ThreadCount = 0;
	double m_MeanQ = 0;
	double m_MeanTC = 0;
	double m_FinalScore = 0;
	bool m_ShowProgress = false;
	vector<double> m_TCs;

public:
	void Load(const string &FileName, const string &RefDir);
	void LoadQ2(const string &FileName, const string &FaDir,
	  const string &RefDirName);
	void AllocThreads();
	double Run(const M3AlnParams &AP);
	void FromSample(const Bench &B, unsigned Pct);
	void Copy(const Bench &B);
	void TCsToFile(const string &FileName) const;
	};
