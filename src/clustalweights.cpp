#include "muscle.h"
#include "clustalweights.h"

void ClustalWeights::Run(const MultiSequence &MS, const Tree &T,
  vector<float> &Weights)
	{
	Weights.clear();
	m_T = &T;
	m_MS = &MS;
	const uint Root = T.GetRoot();
	uint NodeCount = T.GetNodeCount();
	uint SeqCount = MS.GetSeqCount();
	Weights.resize(SeqCount, FLT_MAX);

	m_NodeToSubtreeSize.clear();
	m_NodeToSubtreeSize.resize(NodeCount, 0);

	uint RootSubtreeSize = SetSubtreeSize(Root);
	asserta(RootSubtreeSize == T.GetLeafCount());

	m_NodeToStrength.clear();
	m_NodeToStrength.reserve(NodeCount);
	for (uint Node = 0; Node < NodeCount; ++Node)
		{
		if (T.IsRoot(Node))
			{
			m_NodeToStrength.push_back(0);
			continue;
			}
		uint Parent = T.GetParent(Node);
		float Length = (float) T.GetEdgeLength(Node, Parent);
	// Hack to avoid problems with zero/negative lengths
		if (Length < 0.05f)
			Length = 0.05f;
		uint LeafCount = m_NodeToSubtreeSize[Node];
		float Strength = Length/LeafCount;
		m_NodeToStrength.push_back(Strength);
		}

	float SumWeights = 0;
	for (uint Node = 0; Node < NodeCount; ++Node)
		{
		if (!T.IsLeaf(Node))
			continue;
		uint SeqIndex = T.GetLeafId(Node);
		asserta(SeqIndex < SeqCount && Weights[SeqIndex] == FLT_MAX);
		vector<uint> Path;
		T.GetPathToRoot(Node, Path);
		const uint n = SIZE(Path);
		float Weight = 0;
		for (uint i = 0; i < n; ++i)
			{
			uint Node2 = Path[i];
			float Strength = m_NodeToStrength[Node2];
			Weight += Strength;
			}
		Weights[SeqIndex] = Weight;
		SumWeights += Weight;
		}
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		Weights[SeqIndex] /= SumWeights;
	}

uint ClustalWeights::SetSubtreeSize(uint Node)
	{
	if (m_T->IsLeaf(Node))
		{
		m_NodeToSubtreeSize[Node] = 1;
		return 1;
		}
	uint Left = m_T->GetLeft(Node);
	uint Right = m_T->GetRight(Node);
	uint Size = SetSubtreeSize(Left) + SetSubtreeSize(Right);
	m_NodeToSubtreeSize[Node] = Size;
	return Size;
	}
