#pragma once

#include "tree.h"
#include "pprog.h"
#include "eacluster.h"
#include "treeperm.h"
#include "mpcflat.h"

static const float DEFAULT_MIN_EA_SUPER4_PASS1 = 0.7f;
static const float DEFAULT_MIN_EA_SUPER4_PASS2 = 0.9f;

class Super4
	{
public:
	uint m_TargetPairCount = DEFAULT_TARGET_PAIR_COUNT;
	uint m_MaxClusterSize = DEFAULT_MAX_COARSE_SEQS;
	float m_MinEAPass1 = DEFAULT_MIN_EA_SUPER4_PASS1;
	float m_MinEAPass2 = DEFAULT_MIN_EA_SUPER4_PASS2;

	MultiSequence *m_InputSeqs = 0;

// Pass 1: EA clusters
	EACluster m_EC;
	vector<MultiSequence *> m_ClusterMFAs;
	vector<string> m_ClusterLabels;

// Pass 2: Probcons MSA for each EA cluster
	MPCFlat m_MPC;
	vector<MultiSequence *> m_ClusterMSAs;

// Pass 3: Consensus sequence for each MSA
	MultiSequence m_ConsensusSeqs;

// Pass 4: Distance matrix on consensus sequences
	vector<vector<float> > m_DistMx;

// Pass 5: PProg
	PProg m_PP;
	MultiSequence m_FinalMSA;

	Tree m_GuideTree_None;
	Tree m_GuideTree_ABC;
	Tree m_GuideTree_ACB;
	Tree m_GuideTree_BCA;

	MultiSequence m_FinalMSA_None;
	MultiSequence m_FinalMSA_ABC;
	MultiSequence m_FinalMSA_ACB;
	MultiSequence m_FinalMSA_BCA;

public:
	void SetOpts();
	void ClearTreesAndMSAs();
	void CoarseAlign();
	void InitPP();
	void Run(MultiSequence &InputSeqs, TREEPERM TreePerm);
	uint GetInputSeqCount() const { return m_InputSeqs->GetSeqCount(); }
	void ClusterInput();
	void SplitBigMFA(MultiSequence &MFA, uint MaxSize, float MinEA,
	  vector<MultiSequence *> &SplitMFAs);
	void SplitBigMFA_Random(MultiSequence &MFA, uint MaxSize,
	  vector<MultiSequence *> &SplitMFAs);

	void AlignClusters();
	void GetConsensusSeqs();
	void CalcConsensusSeqsDistMx();
	void MakeGuideTree();
	void DeleteClusterMSAs();
	};
