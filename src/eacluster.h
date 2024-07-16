#pragma once

#include "usorter.h"

class EACluster
	{
public:
	MultiSequence *m_InputSeqs = 0;
	USorter m_US;
	float m_MinEA = FLT_MAX;

	vector<uint> m_CentroidSeqIndexes;
	vector<vector<uint> > m_CentroidIndexToSeqIndexes;
	vector<uint> m_SeqIndexToCentroidIndex;
	vector<MultiSequence *> m_ClusterMFAs;

public:
	void Clear();
	void Run(MultiSequence &InputSeqs, float MinEA);
	void MakeClusterMFAs();
	uint GetBestCentroid(uint SeqIndex, float MinEA, float &BestEA);
	float AlignSeqPair(const string &Label1, const string &Input2);
	void WriteMFAs(const string &FileNamePattern) const;
	void GetClusterMFAs(vector<MultiSequence *> &MFAs) const;
	void Validate() const;
	};
