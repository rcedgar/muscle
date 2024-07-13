#include "muscle.h"
#include "textfile.h"
#include "tree.h"
#include "pprog.h"
#include <map>

void ValidateJoinOrder(const vector<uint> &Indexes1,
  const vector<uint> &Indexes2)
	{
	const uint JoinCount = SIZE(Indexes1);
	asserta(SIZE(Indexes2) == JoinCount);
	
	const uint LeafCount = JoinCount + 1;
	const uint NodeCount = 2*LeafCount - 1;

	set<uint> Pending;
	for (uint LeafIndex = 0; LeafIndex < LeafCount; ++LeafIndex)
		Pending.insert(LeafIndex);

	vector<bool> Used(NodeCount, false);
	for (uint JoinIndex = 0; JoinIndex < JoinCount; ++JoinIndex)
		{
		uint Index1 = Indexes1[JoinIndex];
		uint Index2 = Indexes2[JoinIndex];
		asserta(Index1 != Index2);
		asserta(Index1 < NodeCount);
		asserta(Index2 < NodeCount);

		asserta(Used[Index1] == false);
		asserta(Used[Index2] == false);
		asserta(Pending.find(Index1) != Pending.end());
		asserta(Pending.find(Index2) != Pending.end());

		uint JoinNodeIndex = LeafCount + JoinIndex;
		Used[Index1] = true;
		Used[Index2] = true;

		Pending.erase(Index1);
		Pending.erase(Index2);
		Pending.insert(JoinNodeIndex);
		}
	asserta(SIZE(Pending) == 1);
	uint UsedCount = 0;
	uint NotUsedCount = 0;
	for (uint NodeIndex = 0; NodeIndex < NodeCount; ++NodeIndex)
		{
		if (Used[NodeIndex])
			++UsedCount;
		else
			++NotUsedCount;
		}
	asserta(NotUsedCount == 1);
	}

static const char *GetLabel(const unordered_map<string, uint> &LabelToIndex,
  uint Index)
	{
	for (unordered_map<string, uint>::const_iterator p = LabelToIndex.begin();
	  p != LabelToIndex.end(); ++p)
		{
		if (p->second == Index)
			return p->first.c_str();
		}
	asserta(false);
	return 0;
	}

void LogGuideTreeJoinOrder(const Tree &GuideTree,
  const unordered_map<string, uint> &LabelToIndex,
  const vector<uint> &Indexes1, const vector<uint> &Indexes2)
	{
	asserta(GuideTree.IsRooted());
	const uint NodeCount = GuideTree.GetNodeCount();
	const uint LeafCount = GuideTree.GetLeafCount();
	const uint JoinCount = LeafCount - 1;
	asserta(SIZE(Indexes1) == JoinCount);
	asserta(SIZE(Indexes2) == JoinCount);

	Log("  Join  Index1  Index2\n");
	for (uint JoinIndex = 0; JoinIndex < JoinCount; ++JoinIndex)
		{
		uint Index1 = Indexes1[JoinIndex];
		uint Index2 = Indexes2[JoinIndex];
		Log("%6u", JoinIndex);
		Log("  %6u", Index1);
		Log("  %6u", Index2);

		Log("  ");
		if (Index1 < LeafCount)
			Log(" '%s'", GetLabel(LabelToIndex, Index1));
		else
			Log(" Join%u", Index1 - LeafCount);
		Log(" +");
		if (Index2 < LeafCount)
			Log(" '%s'", GetLabel(LabelToIndex, Index2));
		else
			Log(" Join%u", Index2 - LeafCount);

		Log("\n");
		}
	}

void GetGuideTreeJoinOrder(const Tree &GuideTree,
  const unordered_map<string, uint> &LabelToIndex,
  vector<uint> &Indexes1, vector<uint> &Indexes2)
	{
	asserta(GuideTree.IsRooted());

	Indexes1.clear();
	Indexes2.clear();

	vector<uint> Pending;
	const uint NodeCount = GuideTree.GetNodeCount();
	const uint LeafCount = GuideTree.GetLeafCount();
	vector<bool> IndexUsed(LeafCount);

	const uint JoinCount = LeafCount - 1;
	uint JoinIndex = LeafCount;
	vector<uint> Stack;
	for (uint Node = GuideTree.FirstDepthFirstNode();
	  Node != UINT_MAX;
	  Node = GuideTree.NextDepthFirstNode(Node))
		{
		if (GuideTree.IsLeaf(Node))
			{
			const string Label = GuideTree.GetLeafName(Node);
			unordered_map<string, uint>::const_iterator p = LabelToIndex.find(Label);
			if (p == LabelToIndex.end())
				Die("Label not found >%s", Label.c_str());
			uint Index = p->second;
			asserta(Index < LeafCount);
			asserta(!IndexUsed[Index]);
			Stack.push_back(Index);
			IndexUsed[Index] = true;
			}
		else
			{
			asserta(SIZE(Stack) >= 2);
			uint Left = Stack.back();
			Stack.pop_back();

			uint Right = Stack.back();
			Stack.pop_back();

			Indexes1.push_back(Right);
			Indexes2.push_back(Left);

			Stack.push_back(JoinIndex++);
			}
		}
	}

void MakeGuideTreeFromJoinOrder(
  const vector<uint> &Indexes1, const vector<uint> &Indexes2,
  const unordered_map<string, uint> &LabelToIndex, Tree &GuideTree)
	{
	const uint JoinCount = SIZE(Indexes1);
	asserta(SIZE(Indexes2) == JoinCount);
	const uint LeafCount = JoinCount + 1;
	const uint NodeCount = LeafCount + JoinCount;

	char **LeafLabels = myalloc(char *, LeafCount);
	for (uint LeafIndex = 0; LeafIndex < LeafCount; ++LeafIndex)
		LeafLabels[LeafIndex] = (char *) GetLabel(LabelToIndex, LeafIndex);

	vector<uint> Lefts;
	vector<uint> Rights;
	for (uint JoinIndex = 0; JoinIndex < JoinCount; ++JoinIndex)
		{
		uint Index1 = Indexes1[JoinIndex];
		uint Index2 = Indexes2[JoinIndex];
		Lefts.push_back(Index1);
		Rights.push_back(Index2);
		}

	vector<uint> LeafIds(LeafCount, 1);
	const vector<float> Lengths(NodeCount, 1);

	GuideTree.Create(LeafCount, JoinCount-1,
	  Lefts.data(), Rights.data(), Lengths.data(), Lengths.data(),
	  LeafIds.data(), LeafLabels);
	}

void cmd_guide_tree_join_order()
	{
	const string &TreeFileName = opt(guide_tree_join_order);
	const string &OutputFileName = opt(output);
	Tree GuideTree;
	GuideTree.FromFile(TreeFileName);

	unordered_map<string, uint> LabelToIndex;
	const uint NodeCount = GuideTree.GetNodeCount();
	const uint LeafCount = GuideTree.GetLeafCount();
	const uint JoinCount = LeafCount - 1;

	uint LeafIndex = 0;
	for (uint Node = 0; Node < NodeCount; ++Node)
		{
		if (GuideTree.IsLeaf(Node))
			{
			const string Label = GuideTree.GetLeafName(Node);
			assert(LabelToIndex.find(Label) == LabelToIndex.end());
			LabelToIndex[Label] = LeafIndex++;
			}
		}

	vector<uint> Indexes1;
	vector<uint> Indexes2;
	GetGuideTreeJoinOrder(GuideTree, LabelToIndex, Indexes1, Indexes2);
	LogGuideTreeJoinOrder(GuideTree, LabelToIndex, Indexes1, Indexes2);
	ValidateJoinOrder(Indexes1, Indexes2);

	if (OutputFileName.empty())
		return;

	FILE *f = CreateStdioFile(OutputFileName);
	asserta(GuideTree.IsRooted());
	asserta(SIZE(Indexes1) == JoinCount);
	asserta(SIZE(Indexes2) == JoinCount);

	for (uint JoinIndex = 0; JoinIndex < JoinCount; ++JoinIndex)
		{
		uint Index1 = Indexes1[JoinIndex];
		uint Index2 = Indexes2[JoinIndex];
		fprintf(f, "%u", JoinIndex);
		fprintf(f, "\t%u", Index1);
		fprintf(f, "\t%u", Index2);

		if (Index1 < LeafCount)
			fprintf(f, "\tleaf\t%s", GetLabel(LabelToIndex, Index1));
		else
			fprintf(f, "\tjoin\t%u", Index1 - LeafCount);

		if (Index2 < LeafCount)
			fprintf(f, "\tleaf\t%s", GetLabel(LabelToIndex, Index2));
		else
			fprintf(f, "\tjoin\t%u", Index2 - LeafCount);

		fprintf(f, "\n");
		}

	CloseStdioFile(f);
	}
