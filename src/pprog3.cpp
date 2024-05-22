#include "muscle.h"
#include "pprog3.h"
#include "m3alnparams.h"

void PProg3::Run(const MultiSequence &InputSeqs,
  const vector<float> &InputSeqWeights, const Tree &GuideTree)
	{
	asserta(m_AP != 0);
	const Mx2020 &SubstMx_Letter = m_AP->m_SubstMx_Letter;
	const float GapOpen = m_AP->m_GapOpen;
	m_InputSeqs = &InputSeqs;
	m_InputSeqWeights = &InputSeqWeights;
	m_GuideTree = &GuideTree;
	asserta(GuideTree.IsRooted());
	const uint SeqCount = InputSeqs.GetSeqCount();
	const uint LeafCount = m_GuideTree->GetLeafCount();
	asserta(LeafCount == SeqCount);
	asserta(SIZE(InputSeqWeights) == SeqCount);

	const uint NodeCount = m_GuideTree->GetNodeCount();
	m_NodeToProfile.clear();
	m_NodeToSumInputWeights.clear();
	m_NodeToPath.clear();
	//m_NodeToMSA.clear();

	m_NodeToProfile.resize(NodeCount, 0);
	m_NodeToSumInputWeights.resize(NodeCount, FLT_MAX);
	//m_NodeToMSA.resize(NodeCount, 0);
	m_NodeToPath.resize(NodeCount);

	vector<float> Weights1(1, FLT_MAX);

	uint DoneNodeCount = 0;
	uint Node = GuideTree.FirstDepthFirstNode();
	do
		{
		ProgressStep(DoneNodeCount++, NodeCount, "Progressive align");
		m_NodeToProfile[Node] = 0;
		if (GuideTree.IsLeaf(Node))
			{
			uint SeqIndex = GuideTree.GetLeafId(Node);
			const Sequence &Seq = *InputSeqs.GetSequence(SeqIndex);
#if DEBUG
			{
			string TreeLabel;
			GuideTree.GetLabel(Node, TreeLabel);
			const string &SeqLabel = Seq.GetLabel();
			asserta(TreeLabel == SeqLabel);
			}
#endif
			MultiSequence *MSA1 = new MultiSequence;
			MSA1->AddSequence(&Seq, false);
			float Weight = InputSeqWeights[SeqIndex];
			m_NodeToSumInputWeights[Node] = Weight;

			Weights1[0] = 1.0;
			Profile3 *Prof = new Profile3;
			Prof->FromMSA(*MSA1, SubstMx_Letter, GapOpen, Weights1);
			m_NodeToProfile[Node] = Prof;
			delete MSA1;
			}
		else
			{
			uint Left = GuideTree.GetLeft(Node);
			uint Right = GuideTree.GetRight(Node);
			asserta(Left < NodeCount && Right < NodeCount);

			float SumInputWeightsLeft = m_NodeToSumInputWeights[Left];
			float SumInputWeightsRight = m_NodeToSumInputWeights[Right];
			asserta(SumInputWeightsLeft != FLT_MAX);
			asserta(SumInputWeightsRight != FLT_MAX);
			m_NodeToSumInputWeights[Node] =
			  SumInputWeightsLeft + SumInputWeightsRight;

			Profile3 *ProfLeft = m_NodeToProfile[Left];
			Profile3 *ProfRight = m_NodeToProfile[Right];
			asserta(ProfLeft != 0 && ProfRight != 0);

			string &Path = m_NodeToPath[Node];
			NWSmall3(m_CM, *ProfLeft, *ProfRight, Path);

			if (!GuideTree.IsRoot(Node))
				{
				Profile3 *Prof = new Profile3;
				AlignTwoProfsGivenPath(
				  *ProfLeft, SumInputWeightsLeft,
				  *ProfRight, SumInputWeightsRight,
				  SubstMx_Letter, GapOpen,
				  Path, *Prof);
				m_NodeToProfile[Node] = Prof;
				}

			delete ProfLeft;
			delete ProfRight;
			m_NodeToProfile[Left] = 0;
			m_NodeToProfile[Right] = 0;
			}

		Node = GuideTree.NextDepthFirstNode(Node);
		}
	while (Node != UINT_MAX);
	BuildMSA();
	}

void PProg3::BuildMSA()
	{
	const uint NodeCount = m_GuideTree->GetNodeCount();
	const uint LeafCount = m_GuideTree->GetLeafCount();
	m_MSA.Clear();

	vector<uint> PathToRoot;
	for (uint Node = 0; Node < NodeCount; ++Node)
		{
		ProgressStep(Node, NodeCount, "Root MSA");
		if (!m_GuideTree->IsLeaf(Node))
			continue;
		m_GuideTree->GetPathToRoot(Node, PathToRoot);
		const Sequence *Seq = GetAlignedSeq(Node);
		m_MSA.AddSequence(Seq, true);
		}
	asserta(m_MSA.GetSeqCount() == LeafCount);
	}

const Sequence *PProg3::GetAlignedSeq(uint LeafNode)
	{
	assert(m_GuideTree->IsLeaf(LeafNode));
	uint InputSeqIndex = m_GuideTree->GetLeafId(LeafNode);
	const Sequence *Seq = m_InputSeqs->GetSequence(InputSeqIndex);

	vector<uint> PathToRoot;
	m_GuideTree->GetPathToRoot(LeafNode, PathToRoot);
	const uint N = SIZE(PathToRoot);

// Path includes root, don't got there
	for (uint i = 0; i + 1 < N; ++i)
		{
		uint PathNode = PathToRoot[i];
		uint Parent = m_GuideTree->GetParent(PathNode);
		uint ParentLeft = m_GuideTree->GetLeft(Parent);
		uint ParentRight = m_GuideTree->GetRight(Parent);
		bool IsRight;
		if (ParentRight == PathNode)
			IsRight = true;
		else if (ParentLeft == PathNode)
			IsRight = false;
		else
			asserta(false);
		assert(PathNode < SIZE(m_NodeToPath));
		const string &AlignmentPath = m_NodeToPath[Parent];
		char DorI = (IsRight ? 'I' : 'D');
		const Sequence *ParentSeq = Seq->AddGapsPath(AlignmentPath, DorI);
		if (i != 0)
			DeleteSequence(Seq);
		Seq = ParentSeq;
		}
	return Seq;
	}
