#include "myutils.h"
#include "tree.h"

// ShrubLCAs is non-overlapping subtrees covering entire
// tree such that max size of subtree is n.
void GetShrubs(const Tree &T, uint n, vector<uint> &ShrubLCAs)
	{
	vector<uint> Sizes;
	T.GetSubtreeSizes(Sizes);

	uint ShrubLeafCount = 0;
	const uint NodeCount = T.GetNodeCount();
	const uint LeafCount = T.GetLeafCount();
	for (uint Node = 0; Node < NodeCount; ++Node)
		{
		uint Size = Sizes[Node];
		if (T.IsRoot(Node))
			{
			if (Size <= n)
				{
				asserta(ShrubLCAs.empty());
				ShrubLCAs.push_back(Node);
				ShrubLeafCount = T.GetSubtreeLeafCount(Node);
				break;
				}
			continue;
			}
		uint Parent = T.GetParent(Node);
		uint ParentSize = Sizes[Parent];
		if (Size <= n && ParentSize > n)
			{
			ShrubLCAs.push_back(Node);
			ShrubLeafCount += T.GetSubtreeLeafCount(Node);
			}
		}
	asserta(ShrubLeafCount == LeafCount);
	}

void cmd_shrub()
	{
	Tree T;
	T.FromFile(g_Arg1);
	T.LogMe();
	asserta(T.IsRooted());

	vector<uint> Sizes;
	T.GetSubtreeSizes(Sizes);
	const uint n = (optset_n ? opt_n : 32);

	vector<uint> ShrubLCAs;
	uint ShrubLeafCount = 0;
	const uint NodeCount = T.GetNodeCount();
	const uint LeafCount = T.GetLeafCount();
	for (uint Node = 0; Node < NodeCount; ++Node)
		{
		uint Size = Sizes[Node];
		if (T.IsRoot(Node))
			{
			if (Size <= n)
				{
				asserta(ShrubLCAs.empty());
				ShrubLCAs.push_back(Node);
				ShrubLeafCount = T.GetSubtreeLeafCount(Node);
				break;
				}
			continue;
			}
		uint Parent = T.GetParent(Node);
		uint ParentSize = Sizes[Parent];
		if (Size <= n && ParentSize > n)
			{
			ShrubLCAs.push_back(Node);
			ShrubLeafCount += T.GetSubtreeLeafCount(Node);
			}
		}
	asserta(ShrubLeafCount == LeafCount);
	const uint ShrubCount = SIZE(ShrubLCAs);
	for (uint i = 0; i < ShrubCount; ++i)
		{
		uint LCA = ShrubLCAs[i];
		vector<string> Labels;
		T.GetSubtreeLeafLabels(LCA, Labels);
		const uint M = SIZE(Labels);
		Log("[%4u] %3u ", i, M);
		for (uint j = 0; j < SIZE(Labels); ++j)
			Log(" %s", Labels[j].c_str());
		Log("\n");
		}

	Tree PT;
	vector<string> ShrubLabels;
	PT.PruneTree(T, ShrubLCAs.data(), ShrubCount,
	  "Shrub_", ShrubLabels);

	PT.LogMe();
	uint Node = PT.FirstDepthFirstNode();
	do
		{
		asserta(Node < NodeCount);
		if (!PT.IsLeaf(Node))
			{
			uint Left = PT.GetLeft(Node);
			uint Right = PT.GetRight(Node);
			vector<string> LeftLabels;
			vector<string> RightLabels;
			PT.GetSubtreeLeafLabels(Left, LeftLabels);
			PT.GetSubtreeLeafLabels(Right, RightLabels);
			Log("[");
			bool First = true;
			for (uint i = 0; i < SIZE(LeftLabels); ++i)
				{
				uint L = StrToUint(LeftLabels[i]);
				uint LCA = ShrubLCAs[L];
				vector<string> Labels;
				T.GetSubtreeLeafLabels(LCA, Labels);
				for (uint j = 0; j < SIZE(Labels); ++j)
					{
					if (First)
						First = false;
					else
						Log("+");
					Log("%s", Labels[j].c_str());
					}
				}
			Log("] [");
			for (uint i = 0; i < SIZE(RightLabels); ++i)
				{
				uint L = StrToUint(RightLabels[i]);
				uint LCA = ShrubLCAs[L];
				vector<string> Labels;
				T.GetSubtreeLeafLabels(LCA, Labels);
				for (uint j = 0; j < SIZE(Labels); ++j)
					{
					if (First)
						First = false;
					else
						Log("+");
					Log("%s", Labels[j].c_str());
					}
				}
			Log("]\n");
			}
		Node = PT.NextDepthFirstNode(Node);
		}
	while (Node != UINT_MAX);
	}
