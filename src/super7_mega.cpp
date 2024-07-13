#include "muscle.h"
#include "upgma5.h"
#include "pprog.h"
#include "super7.h"
#include "mpcflat_mega.h"
#include "mega.h"

void CalcGuideTree_SW_BLOSUM62(const MultiSequence &Input, Tree &T);

void Super7_mega::IntraAlignShrub(uint ShrubIndex)
	{
	MPCFlat_mega *MPCm = (MPCFlat_mega *) m_MPC;
	uint LCA = m_ShrubLCAs[ShrubIndex];
	MultiSequence ShrubInput;
	MakeShrubInput(LCA, ShrubInput);

	vector<const vector<vector<byte> > *> ProfilePtrVec;
	GetShrubProfiles(LCA, ProfilePtrVec);

	MPCm->m_TreePerm = TP_None;
	MPCm->Run(&ShrubInput, ProfilePtrVec);
	MultiSequence *ShrubMSA = new MultiSequence;
	ShrubMSA->Copy(*MPCm->m_MSA);
	m_ShrubMSAs.push_back(ShrubMSA);
	}

void Super7_mega::GetShrubProfiles(uint LCA,
  vector<const vector<vector<byte> > *> &ProfilePtrVec)
	{
	ProfilePtrVec.clear();
	vector<uint> LeafNodes;

	m_GuideTree->GetSubtreeLeafNodes(LCA, LeafNodes);
	const uint n = SIZE(LeafNodes);
	asserta(n > 0);
	for (uint i = 0; i < n; ++i)
		{
		uint Node = LeafNodes[i];
		asserta(Node < SIZE(m_NodeToSeqIndex));
		uint SeqIndex = m_NodeToSeqIndex[Node];
		string Label;
		m_GuideTree->GetLabel(Node, Label);
		const vector<vector<byte> > *ptrProfile =
		  Mega::GetProfileByLabel(Label);
		ProfilePtrVec.push_back(ptrProfile);
		}
	}

void cmd_super7_mega()
	{
	uint ShrubSize = 32;
	if (optset_shrub_size)
		ShrubSize = opt(shrub_size);
	if (ShrubSize < 3)
		Die("-shrub_size must be >= 3");

	Mega MM;
	MM.FromFile(g_Arg1);

	MultiSequence InputSeqs;
	InputSeqs.FromStrings(MM.m_Labels, MM.m_Seqs);
	SetGlobalInputMS(InputSeqs);

	bool Nucleo = InputSeqs.GuessIsNucleo();
	CheckMegaOpts(Nucleo);

	Tree GuideTree;
	if (optset_guidetreein)
		GuideTree.FromFile(opt(guidetreein));
	else if (optset_distmxin)
		{
		UPGMA5 U;
		U.ReadDistMx2(opt(distmxin));
		U.Run(LINKAGE_Avg, GuideTree);
		}
	else
		CalcGuideTree_SW_BLOSUM62(InputSeqs, GuideTree);

	SetAlpha(ALPHA_Amino);
	InitProbcons();

	Super7_mega S7;
	MPCFlat_mega *MPCm = new MPCFlat_mega;
	MPCm->m_MM = &MM;
	S7.m_MPC = MPCm;
	S7.Run(InputSeqs, GuideTree, ShrubSize);
	S7.m_FinalMSA.ToFasta(opt(output));

	ProgressLog("Done.\n");
	}
