#pragma once

#include <list>

// ProtDist threshold
class UClustPD
	{
public:
	MultiSequence *m_InputSeqs = 0;
	const vector<uint> *m_SubsetSeqIndexes = 0;
	double m_MaxPD = -1;
	uint m_ThreadCount = 1;
	
	list<uint> m_PendingSubsetIndexes;

// Per-cluster vectors
	vector<uint> m_CentroidSeqIndexes;
	vector<vector<uint> > m_CentroidIndexToMemberSubsetIndexes;

// Per-member vectors
	vector<uint> m_SubsetIndexToCentroidIndex;
	vector<double> m_SubsetIndexToDist;

public:
	void Clear()
		{
		m_ThreadCount = 1;
		m_InputSeqs = 0;
		m_SubsetSeqIndexes = 0;
		m_MaxPD = -1;
		m_PendingSubsetIndexes.clear();
		m_CentroidSeqIndexes.clear();
		m_CentroidIndexToMemberSubsetIndexes.clear();
		m_SubsetIndexToCentroidIndex.clear();
		m_SubsetIndexToDist.clear();
		}

	uint GetInputSeqCount() const { return m_InputSeqs->GetSeqCount(); }
	uint GetSubsetSize() const { return SIZE(*m_SubsetSeqIndexes); }
	const char *GetLabel(uint SeqIndex) const;
	const byte *GetByteSeq(uint SeqIndex, uint &L) const;
	void Run(MultiSequence &InputSeqs,
	  const vector<uint> &SubsetSeqIndexes, double MaxPD);
	double GetProtDistPair(uint SeqIndex1, uint SeqIndex2,
	  string *Path = 0);
	uint Search(uint SeqIndex, const vector<uint> &Centroids,
	  double &BestDist);
	uint SearchAll(uint SeqIndex);
	void ToTsv(FILE *f) const;
	void ToTsv(const string &FileName) const;
	void CentroidsToFasta(const string &FileName) const;
	uint GetClusterCount() const { return SIZE(m_CentroidSeqIndexes); }
	uint GetClusterSize(uint ClusterIndex) const;
	void GetClusterSizes(vector<uint> &Sizes) const;
	void LogStats() const;
	void SelectCandidateGoodCentroids(vector<uint> &SubsetIndexes);
	void GetClusterMFAs(vector<MultiSequence *> &MFAs) const;
	void GetClusterMFA(uint ClusterIndex, MultiSequence &MFA) const;
	};
