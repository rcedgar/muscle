#pragma once

#define DOSC	0

#include "kmerdist66.h"
#include "kmerdist33.h"
#include "clustalweights.h"
#include "pprog3.h"
#include "m3alnparams.h"

#if DOSC
#include "simplecluster.h"
#endif

class Muscle3
	{
public:
	const M3AlnParams *m_AP = 0;
	const MultiSequence *m_InputSeqs = 0;
	vector<vector<float> > m_DistMx;

//K33 AvgQ=0.833 AvgTC=0.533 N=386
//K66 AvgQ=0.835 AvgTC=0.541 N=386
	KmerDist66 m_K66;
	KmerDist33 m_K33;

	vector<string> m_Labels;
	Tree m_GuideTree;
	ClustalWeights m_CW;
	vector<float> m_InputSeqWeights;
	PProg3 m_PP3;
#if DOSC
	SimpleCluster m_SC;
#else
	UPGMA5 m_U5;
#endif
	const MultiSequence *m_FinalMSA = 0;

public:
	void Run(const M3AlnParams &AP, const MultiSequence &InputSeqs);
	void RunRO(const M3AlnParams &AP, const MultiSequence &InputSeqs);
	void WriteMSA(const string &FileName) const;
	};
