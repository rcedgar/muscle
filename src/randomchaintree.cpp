#include "muscle.h"
#include "mpcflat.h"

static void MakeRandomChainTree(const vector<string> &Labels, Tree &T)
	{
	vector<uint> Parents;
	vector<float> Lengths;

	const uint SeqCount = SIZE(Labels);

	vector<uint> SeqIndexes;
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		SeqIndexes.push_back(SeqIndex);
	Shuffle(SeqIndexes);

	Parents.resize(2*SeqCount - 1, UINT_MAX);
	vector<string> NodeLabels;
	for (uint i = 0; i < SeqCount; ++i)
		{
		uint SeqIndex = SeqIndexes[i];
		const string &Label = Labels[SeqIndex];
		Lengths.push_back(1);
		NodeLabels.push_back(Label);
		}

	for (uint i = 0; i + 1 < SeqCount; ++i)
		{
		if (i == 0)
			{
			uint Left = SeqIndexes[0];
			uint Right = SeqIndexes[1];
			Parents[Left] = SeqCount;
			Parents[Right] = SeqCount;
			}
		else
			{
			uint Left = SeqCount + i - 1;
			uint Right = SeqIndexes[i+1];
			Parents[Left] = SeqCount + i;
			Parents[Right] = SeqCount + i;
			}
		NodeLabels.push_back("");
		Lengths.push_back(1);
		}

	T.FromVectors(NodeLabels, Parents, Lengths);
	}

void MPCFlat::CalcGuideTree_RandomChain()
	{
	MakeRandomChainTree(m_Labels, m_GuideTree);
	}

void cmd_labels2randomchaintree()
	{
	const string &LabelsFileName = opt(labels2randomchaintree);
	const string &NewickFileName = opt(output);

	vector<string> Labels;
	ReadStringsFromFile(LabelsFileName, Labels);

	Tree T;
	MakeRandomChainTree(Labels, T);
	T.ToFile(NewickFileName);
	}
