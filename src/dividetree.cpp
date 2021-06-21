#include "muscle.h"
#include "tree.h"

void MakeSubsetNodes(const Tree &InputTree,
  const vector<uint> &SubsetNodes, 
  const vector<string> &SubsetLabels,
  Tree &SubsetTree);

void DivideTree(const Tree &InputTree, uint Node,
  Tree &Subtree, Tree &Supertree)
	{
	asserta(InputTree.IsRooted());
	const uint InputNodeCount = InputTree.GetNodeCount();
	const uint InputLeafCount = InputTree.GetLeafCount();
	asserta(Node < InputNodeCount);
	asserta(!InputTree.IsRoot(Node));

	vector<uint> SubtreeLeafNodes;
	InputTree.GetSubtreeLeafNodes(Node, SubtreeLeafNodes);
	uint N = SIZE(SubtreeLeafNodes);
	asserta(N > 0);

	set<uint> SubtreeSet;
	vector<string> SubtreeLabels;
	for (uint i = 0; i < N; ++i)
		{
		uint Node2 = SubtreeLeafNodes[i];

		string Label;
		InputTree.GetLabel(Node2, Label);

		SubtreeSet.insert(Node2);
		SubtreeLabels.push_back(Label);
		}

	vector<uint> SupertreeLeafNodes;
	vector<string> SupertreeLabels;
	for (uint Node2 = 0; Node2 < InputNodeCount; ++Node2)
		{
		if (!InputTree.IsLeaf(Node2))
			continue;
		if (SubtreeSet.find(Node2) == SubtreeSet.end())
			{
			string Label;
			InputTree.GetLabel(Node2, Label);
			SupertreeLeafNodes.push_back(Node2);
			SupertreeLabels.push_back(Label);
			}
		}

	const uint SubtreeLeafCount = SIZE(SubtreeLeafNodes);
	const uint SupertreeLeafCount = SIZE(SupertreeLeafNodes);
	asserta(SubtreeLeafCount > 0);
	asserta(SupertreeLeafCount > 0);
	asserta(SubtreeLeafCount + SupertreeLeafCount == InputLeafCount);

	MakeSubsetNodes(InputTree, SubtreeLeafNodes, SubtreeLabels,
	  Subtree);

	MakeSubsetNodes(InputTree, SupertreeLeafNodes, SupertreeLabels,
	  Supertree);
	}

void cmd_divide_tree()
	{
	const string &InputFileName = opt(divide_tree);

	Tree InputTree;
	InputTree.FromFile(InputFileName);

	const string &Label1 = opt(label1);
	const string &Label2 = opt(label2);

	uint Node1 = InputTree.GetNodeIndex(Label1);
	uint Node2 = InputTree.GetNodeIndex(Label2);

	uint DivideNode = InputTree.GetLCA(Node1, Node2);

	Tree Subtree;
	Tree Supertree;
	DivideTree(InputTree, DivideNode, Subtree, Supertree);

	Subtree.ToFile(opt(subtreeout));
	Supertree.ToFile(opt(supertreeout));
	}
