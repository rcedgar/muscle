#if 0
#pragma once

#include "tree.h"

class TreeSplitter
	{
public:
	const Tree *m_T = 0;
	uint m_SplitIndex = 0;
	uint m_SplitCount = 0;
	uint m_TargetSize = 0;
	vector<uint> m_SubtreeNodes;

public:
	void Run(const Tree &T, uint SplitCount);
	uint GetBiggestNode() const;
	void GetLabelsVec(vector<vector<string> > &LabelsVec) const;
	void WriteLabels(const string &FileNamePrefix) const;
	void LogState() const;
	void GetSizeOrder(vector<uint> &Order) const;
	void GetSubtree(Tree &Subtree, vector<string> &SplitLabels) const;
	uint GetSplitCount() const { return SIZE(m_SubtreeNodes); }
	};

void MakeSubsetNodes(const Tree &InputTree,
  const vector<uint> &SubsetNodes, 
  const vector<string> &SubsetLabels,
  Tree &SubsetTree);
#endif // 0
