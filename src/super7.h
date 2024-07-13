#pragma once

#include "tree.h"
#include "mpcflat.h"

class Super7
	{
public:
	MultiSequence *m_InputSeqs = 0;
	const Tree *m_GuideTree = 0;
	Tree m_ShrubTree;

	MPCFlat *m_MPC = 0;
	vector<const MultiSequence *> m_ShrubMSAs;
	vector<string> m_ShrubLabels;

	PProg m_PP;
	MultiSequence m_FinalMSA;

	vector<uint> m_NodeToSeqIndex;
	vector<uint> m_SeqIndexToNode;
	vector<uint> m_ShrubLCAs;

public:
	uint GetShrubCount() const { return SIZE(m_ShrubLCAs); }
	void MakeShrubInput(uint LCA, MultiSequence &ShrubInput);

protected:
	void MapLabels();
	void SetShrubs(uint ShrubSize);
	void SetShrubTree();
	void IntraAlignShrubs();
	void ProgAlign();

public:
	virtual void Run(MultiSequence &InputSeqs,
	  const Tree &GuideTree, uint ShrubSize);

private:
	virtual void IntraAlignShrub(uint ShrubIndex);
	};

class Super7_mega : public Super7
	{
public:
	virtual void Run(MultiSequence &InputSeqs,
	  const Tree &GuideTree, uint ShrubSize);

protected:
	void GetShrubProfiles(uint LCA,
	  vector<const vector<vector<byte> > *> &ProfilePtrVec);

protected:
	virtual void IntraAlignShrub(uint ShrubIndex);
	};
