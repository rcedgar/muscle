#include "muscle.h"
#include "upgma5.h"
#include "pprog.h"
#include "super7.h"

void GetShrubs(const Tree &T, uint n, vector<uint> &ShrubLCAs);

void Super7::Run(const MultiSequence &InputSeqs,
  const Tree &GuideTree, uint ShrubSize)
	{
	m_InputSeqs = &InputSeqs;
	m_GuideTree = &GuideTree;
	MapLabels();
	SetShrubs(ShrubSize);
	SetShrubTree();
	IntraAlignShrubs();
	ProgAlign();
	}

void Super7::ProgAlign()
	{
	m_PP.SetMSAs(m_ShrubMSAs, m_ShrubLabels);
	m_PP.RunGuideTree(m_ShrubTree);
	const MultiSequence &FinalMSA = m_PP.GetFinalMSA();
	m_FinalMSA.Copy(FinalMSA);
	//m_FinalMSA.LogMe();
	}

void Super7::SetShrubTree()
	{
	const uint ShrubCount = GetShrubCount();
	m_ShrubTree.PruneTree(*m_GuideTree, m_ShrubLCAs.data(),
	  ShrubCount, "Shrub_", m_ShrubLabels);
	}

void Super7::SetShrubs(uint ShrubSize)
	{
	GetShrubs(*m_GuideTree, ShrubSize, m_ShrubLCAs);
	}

void Super7::MapLabels()
	{
	asserta(m_InputSeqs != 0);
	asserta(m_GuideTree != 0);

	m_SeqIndexToNode.clear();
	m_NodeToSeqIndex.clear();

	const MultiSequence &Seqs = *m_InputSeqs;
	map<string, uint> LabelToSeqIndex;
	const uint SeqCount = Seqs.GetSeqCount();
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const string &Label = Seqs.GetLabelStr(SeqIndex);
		if (LabelToSeqIndex.find(Label) != LabelToSeqIndex.end())
			Die("Duplicate label in sequences >%s", Label.c_str());
		LabelToSeqIndex[Label] = SeqIndex;
		}

	const Tree &T = *m_GuideTree;
	const uint NodeCount = T.GetNodeCount();
	m_NodeToSeqIndex.resize(NodeCount, UINT_MAX);
	m_SeqIndexToNode.resize(SeqCount, UINT_MAX);
	for (uint Node = 0; Node < NodeCount; ++Node)
		{
		if (!T.IsLeaf(Node))
			continue;
		string Label;
		T.GetLabel(Node, Label);
		map<string, uint>::const_iterator iter = LabelToSeqIndex.find(Label);
		if (iter == LabelToSeqIndex.end())
			Die("Tree label not found in sequences >%s", Label.c_str());
		uint SeqIndex = iter->second;
		m_NodeToSeqIndex[Node] = SeqIndex;
		m_SeqIndexToNode[SeqIndex] = Node;
		}
	}

void Super7::MakeShrubInput(uint LCA, MultiSequence &ShrubInput)
	{
	ShrubInput.m_Seqs.clear();
	ShrubInput.m_Owners.clear();

	asserta(m_InputSeqs != 0);
	const MultiSequence &Input = *m_InputSeqs;
	const uint SeqCount = Input.GetSeqCount();
	vector<uint> LeafNodes;
	m_GuideTree->GetSubtreeLeafNodes(LCA, LeafNodes);
	const uint n = SIZE(LeafNodes);
	asserta(n > 0);
	for (uint i = 0; i < n; ++i)
		{
		uint Node = LeafNodes[i];
		asserta(Node < SIZE(m_NodeToSeqIndex));
		uint SeqIndex = m_NodeToSeqIndex[Node];
		asserta(SeqIndex < SeqCount);

		const Sequence *s = Input.m_Seqs[SeqIndex];
		ShrubInput.m_Seqs.push_back(s);
		ShrubInput.m_Owners.push_back(false);
		}
	}

void Super7::IntraAlignShrubs()
	{
	asserta(m_ShrubMSAs.empty());
	uint ShrubCount = GetShrubCount();
	for (uint ShrubIndex = 0; ShrubIndex < ShrubCount; ++ShrubIndex)
		{
		ProgressLog("Aligning shrub %u / %u\n", ShrubIndex+1, ShrubCount);
		uint LCA = m_ShrubLCAs[ShrubIndex];
		MultiSequence ShrubInput;
		MakeShrubInput(LCA, ShrubInput);
		m_MPC.m_TreePerm = TP_None;
		m_MPC.Run(&ShrubInput);
		MultiSequence *ShrubMSA = new MultiSequence;
		ShrubMSA->Copy(*m_MPC.m_MSA);
		m_ShrubMSAs.push_back(ShrubMSA);
		//ShrubMSA->LogMe();
		}
	}

void cmd_super7()
	{
	if (!optset_guidetreein)
		Die("Must set -guidetreein for -super7");

	string &OutputPattern = opt(output);
	if (OutputPattern.empty())
		Die("Must set -output");

	uint ShrubSize = 32;
	if (optset_shrub_size)
		ShrubSize = opt(shrub_size);
	if (ShrubSize < 3)
		Die("-shrub_size must be >= 3");

	Tree GuideTree;
	GuideTree.FromFile(opt(guidetreein));

	LoadGlobalInputMS(g_Arg1);
	ShowGlobalInputSeqStats();

	MultiSequence &InputSeqs = GetGlobalInputMS();
	const uint InputSeqCount = GetGlobalMSSeqCount();

	bool Nucleo = false;
	if (opt(nt))
		Nucleo = true;
	else if (opt(amino))
		Nucleo = false;
	else
		Nucleo = InputSeqs.GuessIsNucleo();

	SetAlpha(Nucleo ? ALPHA_Nucleo : ALPHA_Amino);
	InitProbcons();

	if (optset_diversified)
		Die("-diversified not supported");
	if (optset_replicates)
		Die("-replicates not supported");
	if (optset_stratified)
		Die("-stratified not supported");

	TREEPERM Perm = TP_None;
	if (optset_perm)
		Die("-perm not supported");

	Super7 S7;
	S7.Run(InputSeqs, GuideTree, ShrubSize);
	S7.m_FinalMSA.ToFasta(opt(output));

	ProgressLog("Done.\n");
	}
