#pragma once

class M3AlnParams;

// Progressive alignment from Muscle v3
class PProg3
	{
public:
// Input data
	const M3AlnParams *m_AP = 0;
	const MultiSequence *m_InputSeqs = 0;
	const vector<float> *m_InputSeqWeights = 0;
	const Tree *m_GuideTree = 0;

// Per-node vectors
	vector<Profile3 *> m_NodeToProfile;
	vector<float> m_NodeToSumInputWeights;
	vector<string> m_NodeToPath;
	//vector<MultiSequence *> m_NodeToMSA;
	MultiSequence m_MSA;

// Other
	CacheMem3 m_CM;

public:
	void Run(const MultiSequence &InputSeqs,
	  const vector<float> &InputSeqWeights,
	  const Tree &GuideTree);
	void BuildMSA();
	const Sequence *GetAlignedSeq(uint LeafNode);
	};
