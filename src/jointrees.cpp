#include "muscle.h"
#include "tree.h"

void JoinTrees(const Tree &Tree1, const Tree &Tree2,
  Tree &OutputTree, float NewEdgeLength)
	{
	const uint NodeCount1 = Tree1.GetNodeCount();
	const uint NodeCount2 = Tree2.GetNodeCount();

	vector<string> Labels1;
	vector<string> Labels2;

	vector<uint> Parents1;
	vector<uint> Parents2;

	vector<float> Lengths1;
	vector<float> Lengths2;

	Tree1.ToVectors(Labels1, Parents1, Lengths1);
	Tree2.ToVectors(Labels2, Parents2, Lengths2);

	uint Root = NodeCount1 + NodeCount2;

	vector<string> Labels;
	vector<uint> Parents;
	vector<float> Lengths;

	bool Root1Found = false;
	for (uint Node1 = 0; Node1 < NodeCount1; ++Node1)
		{
		string Label = Labels1[Node1];
		float Length = Lengths1[Node1];
		uint Parent = Parents1[Node1];
		if (Parent == UINT_MAX)
			{
			asserta(!Root1Found);
			Root1Found = true;
			Parent = Root;
			Length = NewEdgeLength;
			}

		Labels.push_back(Label);
		Parents.push_back(Parent);
		Lengths.push_back(Length);
		}
	asserta(Root1Found);

	bool Root2Found = false;
	for (uint Node2 = 0; Node2 < NodeCount2; ++Node2)
		{
		string Label = Labels2[Node2];
		float Length = Lengths2[Node2];
		uint Parent2 = Parents2[Node2];
		uint Parent;
		if (Parent2 == UINT_MAX)
			{
			asserta(!Root2Found);
			Root2Found = true;
			Parent = Root;
			Length = NewEdgeLength;
			}
		else
			Parent = Parent2 + NodeCount1;

		Labels.push_back(Label);
		Parents.push_back(Parent);
		Lengths.push_back(Length);
		}
	asserta(Root2Found);

	Labels.push_back("ROOT");
	Parents.push_back(UINT_MAX);
	Lengths.push_back(0);

	OutputTree.FromVectors(Labels, Parents, Lengths);
	}
