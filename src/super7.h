#pragma once

#include "tree.h"
#include "mpcflat.h"

class Super7
	{
public:
	const MultiSequence *m_InputSeqs = 0;
	const Tree *m_GuideTree = 0;
	Tree m_ShrubTree;

	MPCFlat m_MPC;
	vector<const MultiSequence *> m_ShrubMSAs;
	vector<string> m_ShrubLabels;

	PProg m_PP;
	MultiSequence m_FinalMSA;

	vector<uint> m_NodeToSeqIndex;
	vector<uint> m_SeqIndexToNode;
	vector<uint> m_ShrubLCAs;

public:
	void Run(const MultiSequence &InputSeqs,
	  const Tree &GuideTree, uint ShrubSize);
	uint GetShrubCount() const { return SIZE(m_ShrubLCAs); }
	void MakeShrubInput(uint LCA, MultiSequence &ShrubInput);

private:
	void MapLabels();
	void SetShrubs(uint ShrubSize);
	void SetShrubTree();
	void IntraAlignShrubs();
	void ProgAlign();
	};
