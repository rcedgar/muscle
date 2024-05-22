#pragma once

#include "tree.h"
#include "multisequence.h"

class ClustalWeights
	{
public:
	const Tree *m_T = 0;
	const MultiSequence *m_MS = 0;
	vector<uint> m_NodeToSubtreeSize;
	vector<float> m_NodeToStrength;
	bool m_Dump = false;

public:
	void Run(const MultiSequence &MS, const Tree &T,
	  vector<float> &Weights);
	uint SetSubtreeSize(uint Node);
	};
