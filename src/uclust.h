#pragma once

#include "usorter.h"

class UClust
	{
	const uint MAX_REJECTS = 8;

public:
	MultiSequence *m_InputSeqs = 0;
	float m_MinEA = 0.99f;
	USorter m_US;
	
	vector<uint> m_CentroidSeqIndexes;
	vector<uint> m_SeqIndexToCentroidSeqIndex;
	vector<string> m_SeqIndexToPath;

public:
	void Run(MultiSequence &InputSeqs, float MinEA);
	uint Search(uint SeqIndex, string &Path);
	void AddSeqToIndex(uint SeqIndex);
	float AlignSeqPair(uint SeqIndex1, uint SeqIndex2, string &Path);
	void GetCentroidSeqs(MultiSequence &CentroidSeqs) const;
	void GetGSIs(
	  vector<uint> &CentroidGSIs,
	  vector<uint> &MemberGSIs,
	  vector<uint> &MemberCentroidGSIs,
	  vector<string> &GSIToMemberCentroidPath) const;
	};
