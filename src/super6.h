#pragma once

#include "multisequence.h"
#include "tree.h"
#include "uclustpd.h"
#include "pprog.h"

static const double DEFAULT_MAX_PD_SUPER6_PASS1 = 1.5;
static const uint DEFAULT_TARGET_PAIR_COUNT_CLUSTER_DIST = 8;

class Super6
	{
public:
	double m_MaxPDPass1 = DEFAULT_MAX_PD_SUPER6_PASS1;
	uint m_MaxClusterSize = DEFAULT_MAX_COARSE_SEQS;
	uint m_TargetPairCountClusterDist = DEFAULT_TARGET_PAIR_COUNT_CLUSTER_DIST;
	uint m_TargetPairCount = DEFAULT_TARGET_PAIR_COUNT;
	MultiSequence *m_InputSeqs = 0;

	Tree m_GuideTree;

	MultiSequence m_MSA;
	UClustPD m_UCPD;
	vector<MultiSequence *> m_ClusterMFAs;
	vector<string> m_ClusterLabels;

	MPCFlat m_MPC;
	PProg m_PP;
	vector<MultiSequence *> m_ClusterMSAs;
	vector<vector<float> > m_ClusterDistMx;

public:
	void SetOpts();
	void Run(MultiSequence &InputSeqs);
	void PrepareClusters();
	void AlignClusters();
	//void SplitBigMFA(MultiSequence &MFA, uint MaxSize, float MinEA,
	//  vector<MultiSequence *> &SplitMFAs);
	void SplitBigMFA_Random(MultiSequence &MFA, uint MaxSize,
	  vector<MultiSequence *> &SplitMFAs);
	void CalcClusterDistMx();
	void MakeGuideTree();
	void InitPP();
	};
