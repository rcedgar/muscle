#pragma once

#include "tree.h"
#include "pprog.h"
#include "eacluster.h"
#include "treeperm.h"

class Super4
	{
public:
	MultiSequence *m_InputSeqs = 0;

	string m_SaveDir = "";
	FILE *m_fTSV = 0;

	TREEPERM m_TreePerm = TP_None;

	Tree m_GuideTree;
	//Tree m_GuideTreeABC;
	//Tree m_GuideTreeACB;
	//Tree m_GuideTreeBCA;
	//vector<string> m_LabelsA;
	//vector<string> m_LabelsB;
	//vector<string> m_LabelsC;

// Pass 1: EA clusters
	EACluster m_EC;
	vector<MultiSequence *> m_ClusterMFAs;
	vector<string> m_ClusterLabels;

// Pass 2: Probcons MSA for each EA cluster
	vector<MultiSequence *> m_ClusterMSAs;

// Pass 3: Consensus sequence for each MSA
	MultiSequence m_ConsensusSeqs;

// Pass 4: Distance matrix on consensus sequences
	vector<vector<float> > m_DistMx;

// Pass 5: PProg
	uint m_PairCount = DEFAULT_TARGET_PAIR_COUNT;
	PProg m_PP;
	PProg m_PP_ABC;
	PProg m_PP_ACB;
	PProg m_PP_BCA;
	MultiSequence m_FinalMSA;
	MultiSequence m_FinalMSA_ABC;
	MultiSequence m_FinalMSA_ACB;
	MultiSequence m_FinalMSA_BCA;

public:
	uint GetInputSeqCount() const { return m_InputSeqs->GetSeqCount(); }
	void ClusterInput(uint MaxClusterSize);
	void SplitBigMFA(MultiSequence &MFA, uint MaxSize, float MinEA,
	  vector<MultiSequence *> &SplitMFAs);
	void SplitBigMFA_Random(MultiSequence &MFA, uint MaxSize,
	  vector<MultiSequence *> &SplitMFAs);

	void AlignClusters();
	void GetConsensusSeqs();
	void CalcDistMx();
	void MakeGuideTree();
	void WriteMSAs() const;
	void WriteConsensusSeqs();

	void Final();
	};
