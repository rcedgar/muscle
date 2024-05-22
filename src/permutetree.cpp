#include "muscle.h"
#include "tree.h"
#include "sort.h"

void DivideTree(const Tree &InputTree, uint Node,
  Tree &Subtree, Tree &Supertree);
void JoinTrees(const Tree &Tree1, const Tree &Tree2,
  Tree &OutputTree, float NewEdgeLength);

void StringsToFile(const string &FileName, const vector<string> &v)
	{
	if (FileName.empty())
		return;
	FILE *f = CreateStdioFile(FileName);
	const uint N = SIZE(v);
	for (uint i = 0; i < N; ++i)
		{
		fputs(v[i].c_str(), f);
		fputc('\n', f);
		}
	CloseStdioFile(f);
	}

void DivideTreeFraction(const Tree &InputTree, double Fract,
  Tree &Tree1, Tree &Tree2)
	{
	asserta(Fract > 0 && Fract < 1);

	const uint InputLeafCount = InputTree.GetLeafCount();
	asserta(InputLeafCount >= 3);

	asserta(InputTree.IsRooted());
	const uint NodeCount = InputTree.GetNodeCount();
	uint BestNode = UINT_MAX;
	uint BestLeafCount = UINT_MAX;
	uint BestDiff = UINT_MAX;
	uint TargetLeafCount = uint(InputLeafCount*Fract + 0.5);
	if (TargetLeafCount == 0)
		TargetLeafCount = 1;

	for (uint Node = 0; Node < NodeCount; ++Node)
		{
		uint SubtreeLeafCount = InputTree.GetSubtreeLeafCount(Node);
		uint Diff = (SubtreeLeafCount > TargetLeafCount ?
		  SubtreeLeafCount - TargetLeafCount  :
		  TargetLeafCount  - SubtreeLeafCount);
		if (BestNode == UINT_MAX || Diff < BestDiff)
			{
			BestNode = Node;
			BestDiff = Diff;
			BestLeafCount = SubtreeLeafCount;
			}
		}

	DivideTree(InputTree, BestNode, Tree1, Tree2);
	}

void PermuteTree(const Tree &InputTree,
  Tree &TreeABC, Tree &TreeACB, Tree &TreeBCA,
  vector<string> &LabelsA, vector<string> &LabelsB,
  vector<string> &LabelsC)
	{
	LabelsA.clear();
	LabelsB.clear();
	LabelsC.clear();

	const uint InputLeafCount = InputTree.GetLeafCount();
	if (InputLeafCount < 6)
		{
		TreeABC.Copy(InputTree);
		TreeACB.Copy(InputTree);
		TreeBCA.Copy(InputTree);
		return;
		}

	const float NewEdgeLength = 0.1f;

	Tree TreeA;
	Tree TreeBC;
	Tree TreeB;
	Tree TreeC;
	DivideTreeFraction(InputTree, 0.33, TreeA, TreeBC);
	DivideTreeFraction(TreeBC, 0.5, TreeB, TreeC);

	TreeA.GetLeafLabels(LabelsA);
	TreeB.GetLeafLabels(LabelsB);
	TreeC.GetLeafLabels(LabelsC);

	Tree JoinAB;
	JoinTrees(TreeA, TreeB, JoinAB, NewEdgeLength);
	JoinTrees(JoinAB, TreeC, TreeABC, NewEdgeLength);

	Tree JoinAC;
	JoinTrees(TreeA, TreeC, JoinAC, NewEdgeLength);
	JoinTrees(JoinAC, TreeB, TreeACB, NewEdgeLength);

	Tree JoinBC;
	JoinTrees(TreeB, TreeC, JoinBC, NewEdgeLength);
	JoinTrees(JoinBC, TreeA, TreeBCA, NewEdgeLength);

	TreeABC.Ladderize(true);
	TreeACB.Ladderize(true);
	TreeBCA.Ladderize(true);
	}

void PermTree(Tree &InputTree, TREEPERM TP)
	{
	if (TP == TP_None)
		return;
	uint LeafCount = InputTree.GetLeafCount();
	if (LeafCount < 10)
		return;

	Tree TreeABC;
	Tree TreeACB;
	Tree TreeBCA;
	vector<string> LabelsA;
	vector<string> LabelsB;
	vector<string> LabelsC;
	PermuteTree(InputTree, TreeABC, TreeACB,
	  TreeBCA, LabelsA, LabelsB, LabelsC);
	switch (TP)
		{
	case TP_ABC:
		InputTree.Copy(TreeABC);
		break;

	case TP_ACB:
		InputTree.Copy(TreeACB);
		break;

	case TP_BCA:
		InputTree.Copy(TreeBCA);
		break;

	default:
		asserta(false);
		}
	}

void cmd_permute_tree()
	{
	const string &InputFileName = opt(permute_tree);
	const string &OutputFileName = opt(output);

	Tree InputTree;
	InputTree.FromFile(InputFileName);

	Tree TreeABC;
	Tree TreeACB;
	Tree TreeBCA;
	vector<string> LabelsA;
	vector<string> LabelsB;
	vector<string> LabelsC;
	PermuteTree(InputTree, TreeABC, TreeACB, TreeBCA,
	  LabelsA, LabelsB, LabelsC);

	if (optset_prefix)
		{
		const char *Prefix = opt(prefix).c_str();

		string FileNameABC;
		string FileNameACB;
		string FileNameBCA;

		Ps(FileNameABC, "%sABC.newick", Prefix);
		Ps(FileNameACB, "%sACB.newick", Prefix);
		Ps(FileNameBCA, "%sBCA.newick", Prefix);

		TreeABC.ToFile(FileNameABC);
		TreeACB.ToFile(FileNameACB);
		TreeBCA.ToFile(FileNameBCA);

		StringsToFile(string(Prefix) + "labelsA.txt", LabelsA);
		StringsToFile(string(Prefix) + "labelsB.txt", LabelsB);
		StringsToFile(string(Prefix) + "labelsC.txt", LabelsC);
		}
	}
