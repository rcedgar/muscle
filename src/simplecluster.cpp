#include "myutils.h"
#include "simplecluster.h"

#define TRACE 0

// Everything is indexed by node
// Nodes 0, 1 ... (N-1) are leaves
// Nodes N, N+1 ... (2N-2) are internal nodes
// Total 2N-1 nodes.
// Node (2N-2) is root.

/***
Balibase blosum62
[master dc9b867] Muscle3 with DOSC to use SimpleCluster
-bench d:/a/res/balibase/info/sets.txt -refdir d:/a/res/balibase/ref
  -log bench.log -linkage max

min     AvgQ=0.836 AvgTC=0.537 N=386
biased  AvgQ=0.835 AvgTC=0.538 N=386 b=0.05
biased  AvgQ=0.834 AvgTC=0.535 N=386 b=0.1
biased  AvgQ=0.832 AvgTC=0.534 N=386 b=0.2
avgs    AvgQ=0.820 AvgTC=0.506 N=386
avg     AvgQ=0.818 AvgTC=0.504 N=386
max     AvgQ=0.808 AvgTC=0.480 N=386
***/

float SimpleCluster::GetDist(uint i, uint j) const
	{
	asserta(i < SIZE(m_DistMx));
	asserta(j < SIZE(m_DistMx[i]));
	float d = m_DistMx[i][j];
	asserta(feq(d, m_DistMx[j][i]));
	return d;
	}

float SimpleCluster::FindClosestPair(uint &Index1, uint &Index2) const
	{
	Index1 = UINT_MAX;
	Index2 = UINT_MAX;
	asserta(!m_Pending.empty());
	float BestDist = FLT_MAX;
	for (list<uint>::const_iterator iter_i = m_Pending.begin();
	  iter_i != m_Pending.end(); ++iter_i)
		{
		uint i = *iter_i;
		list<uint>::const_iterator iter_j = iter_i;
		++iter_j;
		while (iter_j != m_Pending.end())
			{
			uint j = *iter_j;
			float d = GetDist(i, j);
			if (BestDist == FLT_MAX)
				{
				Index1 = i;
				Index2 = j;
				BestDist = d;
				}
			else
				{
				bool Closer = (m_DistIsSimilarity ?
				  (d > BestDist) : (d < BestDist));
				if (Closer)
					{
					Index1 = i;
					Index2 = j;
					BestDist = d;
					}
				}
			++iter_j;
			}
		}
	asserta(Index1 != UINT_MAX && Index2 != UINT_MAX);
	return BestDist;
	}

uint SimpleCluster::GetSize(uint i) const
	{
	asserta(i < SIZE(m_Sizes));
	uint n = m_Sizes[i];
	return n;
	}

// i1,i2 are nearest neighbors forming new cluster
// j is existing cluster
float SimpleCluster::CalcNewDist(uint i1, uint i2, uint j) const
	{
	asserta(i1 < SIZE(m_DistMx));
	asserta(i2 < SIZE(m_DistMx));
	asserta(j < SIZE(m_DistMx));
	asserta(i1 != i2 && i1 != j && i2 != j);

	uint Size1 = GetSize(i1);
	uint Size2 = GetSize(i2);
	asserta(Size1 > 0 && Size2 > 0);

	float d1 = GetDist(i1, j);
	float d2 = GetDist(i2, j);

	float NewDist = FLT_MAX;
	if (m_Linkage == "avg")
		NewDist = (d1 + d2)/2;
	else if (m_Linkage == "avgs")
		NewDist = (d1*Size1 + d2*Size2)/(Size1 + Size2);
	else if (m_Linkage == "min")
		NewDist = min(d1,  d2);
	else if (m_Linkage == "max")
		NewDist = max(d1,  d2);
	else if (m_Linkage == "biased")
		{
		float b = 0.05f;
		NewDist = b*(d1 + d2)/2 + (1 - b)*min(d1, d2);
		}
	else
		Die("SimpleCluster::m_Linkage=%s", m_Linkage.c_str());

	return NewDist;
	}

void SimpleCluster::Join(uint JoinIndex)
	{
	asserta(m_Pending.size() > 1);
	
	uint Node = m_InputCount + JoinIndex;

	uint Index1, Index2;
	FindClosestPair(Index1, Index2);

	list<uint>::const_iterator iterIndex1 = m_Pending.end();
	list<uint>::const_iterator iterIndex2 = m_Pending.end();
	for (list<uint>::const_iterator iter = m_Pending.begin();
	  iter != m_Pending.end(); ++iter)
		{
		uint Index3 = *iter;
		if (Index3 == Index1)
			{
			iterIndex1 = iter;
			continue;
			}
		else if (Index3 == Index2)
			{
			iterIndex2 = iter;
			continue;
			}
		float NewDist = CalcNewDist(Index1, Index2, Index3);
		asserta(m_DistMx[Node][Index3] == FLT_MAX);
		asserta(m_DistMx[Index3][Node] == FLT_MAX);
		m_DistMx[Node][Index3] = NewDist;
		m_DistMx[Index3][Node] = NewDist;
		}

	uint Size1 = GetSize(Index1);
	uint Size2 = GetSize(Index2);
	asserta(m_Sizes[Node] == 0);
	m_Sizes[Node] = Size1 + Size2;

	m_Parents[Index1] = Node;
	m_Parents[Index2] = Node;

	m_Lefts[Node] = Index1;
	m_Rights[Node] = Index2;

	float d12 = m_DistMx[Index1][Index2];
	float Height = d12/2;
	m_Heights[Node] = Height;

	float Height1 = m_Heights[Index1];
	float Height2 = m_Heights[Index2];

	m_Lengths[Index1] = Height - Height1;
	m_Lengths[Index2] = Height - Height2;
	
	m_Pending.erase(iterIndex1);
	m_Pending.erase(iterIndex2);
	m_Pending.push_back(Node);
	}

void SimpleCluster::Run(const vector<vector<float> > &DistMx,
  const vector<string> &Labels, const string &Linkage,
  bool DistIsSimilarity)
	{
	Clear();
	m_DistMx = DistMx;
	m_Labels = Labels;
	m_Linkage = Linkage;
	m_DistIsSimilarity = DistIsSimilarity;

	m_InputCount = SIZE(Labels);
	const uint NodeCount = 2*m_InputCount - 1;
	const uint JoinCount = m_InputCount - 1;

	m_DistMx.resize(NodeCount);
	for (uint Node = 0; Node < NodeCount; ++Node)
		m_DistMx[Node].resize(NodeCount, FLT_MAX);

	for (uint JoinIndex = 0; JoinIndex < JoinCount; ++JoinIndex)
		{
		string Label;
		Ps(Label, "Int%u", JoinIndex);
		m_Labels.push_back(Label);
		}

	m_Parents.resize(NodeCount, UINT_MAX);
	m_Lefts.resize(NodeCount, UINT_MAX);
	m_Rights.resize(NodeCount, UINT_MAX);
	m_Sizes.resize(NodeCount, 0);
	m_Lengths.resize(NodeCount, FLT_MAX);
	m_Heights.resize(NodeCount, FLT_MAX);

	for (uint Node = 0; Node < m_InputCount; ++Node)
		{
		m_Pending.push_back(Node);
		m_Sizes[Node] = 1;
		m_Heights[Node] = 0;
		}

	for (uint JoinIndex = 0; JoinIndex < JoinCount; ++JoinIndex)
		{
#if TRACE
		Log("\n\n___________ Join %u ___________\n", JoinIndex);
		LogMe();
#endif
		Join(JoinIndex);
		}
#if TRACE
	LogMe();
#endif
	}

const string &SimpleCluster::GetLabel(uint Node) const
	{
	asserta(Node < SIZE(m_Labels));
	return m_Labels[Node];
	}

void SimpleCluster::LogMe() const
	{
	Log("%8.8s", "DistMx");
	const uint NodeCount = SIZE(m_DistMx);
	for (uint Node = 0; Node < NodeCount; ++Node)
		Log("  %8.8s", GetLabel(Node).c_str());
	Log("\n");

	for (uint Node = 0; Node < NodeCount; ++Node)
		{
		const string &Label = GetLabel(Node);
		Log("%8.8s", Label.c_str());
		for (uint Node2 = 0; Node2 < NodeCount; ++Node2)
			{
			if (Node == Node2)
				Log("  %8.8s", ".");
			else
				{
				float d = GetDist(Node, Node2);
				if (d == FLT_MAX)
					Log("  %8.8s", "*");
				else
					Log("  %8.3g", d);
				}
			}
		Log("\n");
		}
	Log("Pending (%u) ", SIZE(m_Pending));
	for (list<uint>::const_iterator iter = m_Pending.begin();
	  iter != m_Pending.end(); ++iter)
		Log(" %u", *iter);
	Log("\n");
	Log("  Node    Size  Height  Length  Parent    Left   Right  Label\n");
	//   123456  123456  123456  123456  123456  123456  123456  123456
	for (uint Node = 0; Node < NodeCount; ++Node)
		{
		uint Parent = m_Parents[Node];
		uint Left = m_Lefts[Node];
		uint Right = m_Rights[Node];
		float Height = m_Heights[Node];
		float Length = m_Lengths[Node];
		const string &Label = GetLabel(Node);
		Log("%6u", Node);
		Log("  %6u", m_Sizes[Node]);
		if (Height == FLT_MAX) Log("  %6.6s", "*"); else Log("  %6.3g", Height);
		if (Length == FLT_MAX) Log("  %6.6s", "*"); else Log("  %6.3g", Length);
		if (Parent == UINT_MAX) Log("  %6.6s", "*"); else Log("  %6u", Parent);
		if (Left == UINT_MAX) Log("  %6.6s", "*"); else Log("  %6u", Left);
		if (Right == UINT_MAX) Log("  %6.6s", "*"); else Log("  %6u", Right);
		Log("  >%s", Label.c_str());
		Log("\n");
		}
	}

void SimpleCluster::GetTree(Tree &T) const
	{
	T.FromVectors(m_Labels, m_Parents, m_Lengths);
	}
