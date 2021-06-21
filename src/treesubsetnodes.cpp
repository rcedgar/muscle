#include "myutils.h"
#include "muscle.h"
#include "tree.h"
#include <set>

void IntsFromFile(const string &FileName, vector<uint> &Ints,
  vector<string> &Labels)
	{
	Ints.clear();
	Labels.clear();
	FILE *f = OpenStdioFile(FileName);
	string Line;
	vector<string> Fields;
	Progress("Reading %s...", FileName.c_str());
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 2);
		uint i = StrToUint(Fields[0]);
		const string &Label = Fields[1];
		Ints.push_back(i);
		Labels.push_back(Label);
		}
	Progress("done.\n");
	CloseStdioFile(f);
	}

void MakeSubsetNodes(const Tree &InputTree,
  const vector<uint> &SubsetNodes, 
  const vector<string> &SubsetLabels,
  Tree &SubsetTree)
	{
	if (!InputTree.IsRooted())
		Die("Tree must be rooted");

	const uint SubsetNodeCount = SIZE(SubsetNodes);
	if (SubsetNodeCount < 2)
		Die("Need at least two nodes");
	asserta(SIZE(SubsetLabels) == SubsetNodeCount);

	const uint NodeCount = InputTree.GetNodeCount();

	map<uint, string> OldNodeToNewLabel;
	set<uint> SubsetSet;
	for (uint i = 0; i < SubsetNodeCount; ++i)
		{
		uint SubsetNode = SubsetNodes[i];
		SubsetSet.insert(SubsetNode);
		OldNodeToNewLabel[SubsetNode] = SubsetLabels[i];
		}

// ParentSet includes subset nodes which are
// parents of other subset nodes
	set<uint> ParentSet;
	vector<uint> Path;
	for (uint i = 0; i < SubsetNodeCount; ++i)
		{
		uint SubsetNode = SubsetNodes[i];
		InputTree.GetPathToRoot(SubsetNode, Path);
		const uint n = SIZE(Path);
		asserta(Path[0] == SubsetNode);
		for (uint j = 1; j < n; ++j)
			{
			uint Node = Path[j];
			asserta(Node < NodeCount);
			if (SubsetSet.find(Node) != SubsetSet.end())
				{
				ProgressLog("Node %u is parent of %u\n",
				  Node, SubsetNode);
				ParentSet.insert(Node);
				}
			}
		}
	const uint ParentCount = SIZE(ParentSet);
	set<uint> NewTips;
	for (set<uint>::const_iterator p = SubsetSet.begin();
	  p != SubsetSet.end(); ++p)
		{
		uint Node = *p;
		if (ParentSet.find(Node) == ParentSet.end())
			NewTips.insert(Node);
		}
	const uint NewTipCount = SIZE(NewTips);
	if (NewTipCount == 1)
		Die("One tip in subset tree");

	asserta(NewTipCount + ParentCount == SubsetNodeCount);

	vector<uint> NodeToPathCount(NodeCount);
	for (set<uint>::const_iterator p = NewTips.begin(); p != NewTips.end(); ++p)
		{
		uint Tip = *p;
		InputTree.GetPathToRoot(Tip, Path);
		const uint n = SIZE(Path);
		asserta(Path[0] == Tip);
		for (uint j = 0; j < n; ++j)
			{
			uint Node2 = Path[j];
			++NodeToPathCount[Node2];
			}
		}

	vector<uint> OldNodeToNewParent(NodeCount, UINT_MAX);
	set<uint> Pending;
	set<uint> Done;
	for (set<uint>::const_iterator p = NewTips.begin();
	  p != NewTips.end(); ++p)
		Pending.insert(*p);

	for (;;)
		{
		if (Pending.empty())
			break;
		uint Node = *Pending.begin();
		asserta(Done.find(Node) == Done.end());
		Done.insert(Node);
		Pending.erase(Node);
		if (InputTree.IsRoot(Node))
			continue;

		InputTree.GetPathToRoot(Node, Path);
		asserta(Path[0] == Node);
		const uint n = SIZE(Path);
		uint NewParent = UINT_MAX;
		for (uint j = 1; j < n; ++j)
			{
			uint Node2 = Path[j];
			uint Left = InputTree.GetLeft(Node2);
			uint Right = InputTree.GetRight(Node2);
			asserta(Left != UINT_MAX);
			asserta(Right != UINT_MAX);
			uint LeftPathCount = NodeToPathCount[Left];
			uint RightPathCount = NodeToPathCount[Right];
			if (LeftPathCount > 0 && RightPathCount > 0)
				{
				NewParent = Node2;
				break;
				}
			}
		asserta(OldNodeToNewParent[Node] == UINT_MAX);
		OldNodeToNewParent[Node] = NewParent;
		if (NewParent != UINT_MAX && Done.find(NewParent) == Done.end())
			Pending.insert(NewParent);
		}
	
	vector<uint> OldNodeToNewNode(NodeCount, UINT_MAX);
	vector<uint> NewNodeToOldNode;
	for (set<uint>::const_iterator p = Done.begin(); p != Done.end(); ++p)
		{
		uint OldNode = *p;
		asserta(OldNode < SIZE(OldNodeToNewNode));
		uint NewNode = SIZE(NewNodeToOldNode);

		OldNodeToNewNode[OldNode] = NewNode;
		NewNodeToOldNode.push_back(OldNode);
		}

	vector<uint> NewParents;
	vector<float> NewLengths;
	vector<string> NewLabels;
	for (set<uint>::const_iterator p = Done.begin(); p != Done.end(); ++p)
		{
		uint OldNode = *p;
		asserta(OldNode < SIZE(OldNodeToNewNode));
		uint NewNode = SIZE(NewNodeToOldNode);
		asserta(OldNode < SIZE(OldNodeToNewParent));
		uint Parent = OldNodeToNewParent[OldNode];
		uint NewParent = UINT_MAX;
		if (Parent != UINT_MAX)
			NewParent = OldNodeToNewNode[Parent];

		string Label;
		map<uint, string>::const_iterator q = OldNodeToNewLabel.find(OldNode);
		if (q == OldNodeToNewLabel.end())
			InputTree.GetLabel(OldNode, Label);
		else
			Label = q->second;

		float Distance = (Parent == UINT_MAX ? 0 : (float) InputTree.GetDistance(OldNode, Parent));

		NewParents.push_back(NewParent);
		NewLabels.push_back(Label);
		NewLengths.push_back(Distance);
		}

	SubsetTree.FromVectors(NewLabels, NewParents, NewLengths);
	}

void cmd_tree_subset_nodes()
	{
	const string &InputFileName = opt(tree_subset_nodes);
	const string &NodesFileName = opt(nodes);

	vector<uint> Nodes;
	vector<string> NewLabels;
	IntsFromFile(NodesFileName, Nodes, NewLabels);
	
	Tree T;
	T.FromFile(InputFileName);

	Tree Subtree;
	MakeSubsetNodes(T, Nodes, NewLabels, Subtree);
	T.Ladderize(opt(right));
	Subtree.ToFile(opt(output));
	}
