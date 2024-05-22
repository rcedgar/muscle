#pragma once

#include <list>
#include "tree.h"

class SimpleCluster
	{
public:
	uint m_InputCount = 0;
	vector<vector<float> > m_DistMx;
	vector<string> m_Labels;
	vector<uint> m_Sizes;
	string m_Linkage;
	bool m_DistIsSimilarity = false;

	list<uint> m_Pending;
	vector<uint> m_Parents;
	vector<uint> m_Lefts;
	vector<uint> m_Rights;
	vector<float> m_Lengths;
	vector<float> m_Heights;

public:
	void Clear()
		{
		m_InputCount = 0;
		m_DistMx.clear();
		m_Labels.clear();
		m_Sizes.clear();
		m_Linkage = "Undefined";
		m_DistIsSimilarity = false;
		m_Pending.clear();
		m_Parents.clear();
		m_Lefts.clear();
		m_Rights.clear();
		m_Lengths.clear();
		m_Heights.clear();
		}

	void Run(const vector<vector<float> > &DistMx,
	  const vector<string> &Labels, const string &Linkage,
	  bool DistIsSimilarity);
	const string &GetLabel(uint Node) const;
	float FindClosestPair(uint &Index1, uint &Index2) const;
	float GetDist(uint i, uint j) const;
	float CalcNewDist(uint i1, uint i2, uint j) const;
	uint GetSize(uint i) const;
	void Join(uint JoinIndex);
	void LogMe() const;
	void GetTree(Tree &T) const;
	};
