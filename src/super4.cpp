#include "muscle.h"
#include "super4.h"
#include "sequence.h"
#include "multisequence.h"
#include "usorter.h"
#include "upgma5.h"
#include "pprog.h"
#include "treeperm.h"

void LogDistMx(const string &Msg, const vector<vector<float> > &Mx);
void GetConsensusSequence(const MultiSequence &MSA, string &Seq);

void Super4::ClearTreesAndMSAs()
	{
	m_FinalMSA.Clear();

	m_GuideTree_None.Clear();
	m_GuideTree_ABC.Clear();
	m_GuideTree_ACB.Clear();
	m_GuideTree_BCA.Clear();

	m_FinalMSA_None.Clear();
	m_FinalMSA_ABC.Clear();
	m_FinalMSA_ACB.Clear();
	m_FinalMSA_BCA.Clear();
	}

void Super4::MakeGuideTree()
	{
	UPGMA5 U;
	U.Init(m_ClusterLabels, m_DistMx);
	U.FixEADistMx();
	U.Run(LINKAGE_Biased, m_GuideTree_None);
	}

void Super4::SplitBigMFA_Random(MultiSequence &InputMFA, uint MaxSize,
  vector<MultiSequence *> &SplitMFAs)
	{
	SplitMFAs.clear();
	const uint InputSeqCount = InputMFA.GetSeqCount();
	asserta(InputSeqCount > MaxSize);
	uint OutputSeqCount = 0;
	for (;;)
		{
		asserta(OutputSeqCount <= InputSeqCount);
		uint RemainingSeqCount = InputSeqCount - OutputSeqCount;
		if (RemainingSeqCount == 0)
			break;
		
		uint N = RemainingSeqCount;
		if (N > MaxSize)
			N = MaxSize;

		MultiSequence *SplitMFA = new MultiSequence;
		asserta(SplitMFA != 0);
		for (uint i = 0; i < N; ++i)
			{
			const Sequence *seq = InputMFA.GetSequence(OutputSeqCount + i);
			SplitMFA->AddSequence(seq, false);
			}
		SplitMFAs.push_back(SplitMFA);

		OutputSeqCount += N;
		}
	AssertSameSeqsVec(InputMFA, SplitMFAs);
	}

void Super4::SplitBigMFA(MultiSequence &BigMFA, uint MaxSize, float MinEA,
  vector<MultiSequence *> &SplitMFAs)
	{
	SplitMFAs.clear();

	const uint InputSeqCount = BigMFA.GetSeqCount();
	asserta(InputSeqCount > MaxSize);

	m_EC.Run(BigMFA, MinEA);
	m_EC.GetClusterMFAs(SplitMFAs);
	uint ClusterCount = SIZE(SplitMFAs);
	asserta(ClusterCount > 0);
	AssertSameSeqsVec(BigMFA, SplitMFAs);

	for (uint ClusterIndex = 0; ClusterIndex < ClusterCount; ++ClusterIndex)
		{
		MultiSequence *MFA = SplitMFAs[ClusterIndex];
		uint SeqCount = MFA->GetSeqCount();
		if (SeqCount > MaxSize)
			{
			vector<MultiSequence *> *SubMFAs = new vector<MultiSequence *>;
			SplitBigMFA_Random(*MFA, MaxSize, *SubMFAs);
			AssertSameSeqsVec(*MFA, *SubMFAs);

			const uint N = SIZE(*SubMFAs);
			asserta(N > 1);
			SplitMFAs[ClusterIndex] = (*SubMFAs)[0];
			for (uint i = 1; i < N; ++i)
				{
				MultiSequence *SubMFA = (*SubMFAs)[i];
				SplitMFAs.push_back(SubMFA);
				}
			}
		}
	AssertSameSeqsVec(BigMFA, SplitMFAs);
	}

void Super4::ClusterInput()
	{
	const uint InputSeqCount = m_InputSeqs->GetSeqCount();

	m_EC.Run(*m_InputSeqs, m_MinEAPass1);
	m_EC.GetClusterMFAs(m_ClusterMFAs);
	AssertSameSeqsVec(*m_InputSeqs, m_ClusterMFAs);

	uint ClusterCount = SIZE(m_ClusterMFAs);
	ProgressLog("%u clusters pass 1\n", ClusterCount);
	for (uint ClusterIndex = 0; ClusterIndex < ClusterCount; ++ClusterIndex)
		{
		MultiSequence *MFA = m_ClusterMFAs[ClusterIndex];
		AssertSameLabels(*MFA);
		uint SeqCount = MFA->GetSeqCount();
		asserta(SeqCount > 0);
		if (SeqCount > m_MaxClusterSize)
			{
			vector<MultiSequence *> *SplitMFAs = new vector<MultiSequence *>;
			SplitBigMFA(*MFA, m_MaxClusterSize, m_MinEAPass2, *SplitMFAs);

			const uint N = SIZE(*SplitMFAs);
			asserta(N > 1);

			m_ClusterMFAs[ClusterIndex] = (*SplitMFAs)[0];
			for (uint i = 1; i < N; ++i)
				{
				MultiSequence *SplitMFA = (*SplitMFAs)[i];
				m_ClusterMFAs.push_back(SplitMFA);
				AssertSameLabels(*SplitMFA);
				}
			}
		}

	ClusterCount = SIZE(m_ClusterMFAs);
	ProgressLog("%u clusters pass 2\n", ClusterCount);

	AssertSameSeqsVec(*m_InputSeqs, m_ClusterMFAs);

	m_ClusterLabels.clear();
	for (uint ClusterIndex = 0; ClusterIndex < ClusterCount; ++ClusterIndex)
		{
		MultiSequence *MFA = m_ClusterMFAs[ClusterIndex];
		uint SeqCount = MFA->GetSeqCount();
		asserta(SeqCount <= m_MaxClusterSize);

		string ClusterLabel;
		Ps(ClusterLabel, "Cluster%u", ClusterIndex);
		m_ClusterLabels.push_back(ClusterLabel);
		}
	}

void Super4::AlignClusters()
	{
	m_ClusterMSAs.clear();

	uint ClusterCount = SIZE(m_ClusterMFAs);
	for (uint ClusterIndex = 0; ClusterIndex < ClusterCount; ++ClusterIndex)
		{
		MultiSequence *ClusterMFA = m_ClusterMFAs[ClusterIndex];
		AssertSameLabels(*ClusterMFA);
		const uint SeqCount = ClusterMFA->GetSeqCount();

		if (SeqCount == 1)
			Progress("Align cluster %u / %u (1 seq)\n",
			  ClusterIndex + 1, ClusterCount);
		else
			{
			Progress("\n");
			Progress("Align cluster %u / %u (%u seqs)\n",
			  ClusterIndex + 1, ClusterCount, SeqCount);
			Progress("\n");
			}

		m_MPC.m_TreePerm = TP_None;
		m_MPC.Run(ClusterMFA);
		MultiSequence *ClusterMSA = new MultiSequence;
		asserta(ClusterMSA != 0);
		ClusterMSA->Copy(*m_MPC.m_MSA);
		AssertSameLabels(*ClusterMSA);
		m_ClusterMSAs.push_back(ClusterMSA);
		}
	}

void Super4::DeleteClusterMSAs()
	{
	const uint N = SIZE(m_ClusterMSAs);
	for (uint i = 0; i < N; ++i)
		{
		MultiSequence *MSA = m_ClusterMSAs[i];
		if (MSA != 0)
			MSA->Clear();
		}
	m_ClusterMSAs.clear();
	}

void Super4::GetConsensusSeqs()
	{
	m_ConsensusSeqs.Clear();
	uint ClusterCount = SIZE(m_ClusterMSAs);
	for (uint ClusterIndex = 0; ClusterIndex < ClusterCount; ++ClusterIndex)
		{
		ProgressStep(ClusterIndex, ClusterCount, "Consensus sequences");
		MultiSequence *ClusterMSA = m_ClusterMSAs[ClusterIndex];

		string Label;
		Ps(Label, "Cluster%u", ClusterIndex);

		string ConsSeq;
		GetConsensusSequence(*ClusterMSA, ConsSeq);

		Sequence *seq = NewSequence();
		seq->FromString(Label, ConsSeq);
		
		AddGlobalTmpSeq(seq);
		m_ConsensusSeqs.AddSequence(seq, true);

		asserta(ClusterIndex < SIZE(m_ClusterLabels));
		const string &ClusterLabel = m_ClusterLabels[ClusterIndex];
		}
	if (optset_calnout)
		{
		m_MPC.Run(&m_ConsensusSeqs);
		m_MPC.m_MSA->WriteMFA(opt(calnout));
		}
	}

void Super4::InitPP()
	{
	vector<const MultiSequence *> MSAs;
	const uint n = SIZE(m_ClusterMSAs);
	for (uint i = 0; i < n; ++i)
		{
		const MultiSequence *MSA = m_ClusterMSAs[i];
		MSAs.push_back(MSA);
		}

	m_PP.m_TargetPairCount = m_TargetPairCount;
	m_PP.SetMSAs(MSAs, m_ClusterLabels);
	}

void Super4::SetOpts()
	{
	m_TargetPairCount = optd(paircount, DEFAULT_TARGET_PAIR_COUNT);
	m_MaxClusterSize = optd(paircount, DEFAULT_MAX_COARSE_SEQS);
	m_MinEAPass1 = (float) optd(super4_minea1, DEFAULT_MIN_EA_SUPER4_PASS1);
	m_MinEAPass2 = (float) optd(super4_minea2, DEFAULT_MIN_EA_SUPER4_PASS2);
	m_MPC.m_ConsistencyIterCount = optd(consiters, 2);
	m_MPC.m_RefineIterCount = optd(refineiters, 100);
	}

void Super4::CalcConsensusSeqsDistMx()
	{
	FILE *f = 0;
	CalcEADistMx(f, &m_ConsensusSeqs, m_DistMx);
	}

void Super4::CoarseAlign()
	{
	AssertSameLabels(*m_InputSeqs);
	ClusterInput();
	AlignClusters();
	GetConsensusSeqs();
	CalcConsensusSeqsDistMx();
	MakeGuideTree();
	InitPP();
	}

void Super4::Run(MultiSequence &InputSeqs, TREEPERM TreePerm)
	{
	m_InputSeqs = &InputSeqs;
	SetOpts();
	CoarseAlign();

	if (TreePerm == TP_None)
		{
		m_PP.RunGuideTree(m_GuideTree_None);
		const MultiSequence &FinalMSA = m_PP.GetFinalMSA();
		m_FinalMSA.Copy(FinalMSA);
		DeleteClusterMSAs();
		return;
		}

	vector<string> LabelsA;
	vector<string> LabelsB;
	vector<string> LabelsC;
	PermuteTree(m_GuideTree_None,
	  m_GuideTree_ABC, m_GuideTree_ACB, m_GuideTree_BCA,
	  LabelsA, LabelsB, LabelsC);

	switch (TreePerm)
		{
	case TP_ABC:
		ProgressLog("Guide tree ABC\n");
		m_PP.RunGuideTree(m_GuideTree_ABC);
		m_FinalMSA.Copy(m_PP.GetFinalMSA());
		break;

	case TP_ACB:
		ProgressLog("Guide tree ACB\n");
		m_PP.RunGuideTree(m_GuideTree_ACB);
		m_FinalMSA.Copy(m_PP.GetFinalMSA());
		break;

	case TP_BCA:
		ProgressLog("Guide tree BCA\n");
		m_PP.RunGuideTree(m_GuideTree_BCA);
		m_FinalMSA.Copy(m_PP.GetFinalMSA());
		break;

	case TP_All:
		ProgressLog("Guide tree (default)\n");
		m_PP.RunGuideTree(m_GuideTree_None);
		m_FinalMSA_None.Copy(m_PP.GetFinalMSA());

		ProgressLog("Guide tree ABC\n");
		m_PP.RunGuideTree(m_GuideTree_ABC);
		m_FinalMSA_ABC.Copy(m_PP.GetFinalMSA());

		ProgressLog("Guide tree ACB\n");
		m_PP.RunGuideTree(m_GuideTree_ACB);
		m_FinalMSA_ACB.Copy(m_PP.GetFinalMSA());

		ProgressLog("Guide tree BCA\n");
		m_PP.RunGuideTree(m_GuideTree_BCA);
		m_FinalMSA_BCA.Copy(m_PP.GetFinalMSA());
		break;

	default:
		asserta(false);
		}
	DeleteClusterMSAs();
	}

void cmd_super4()
	{
	if (optset_mega)
		Die("-super4 does not support -mega, use -super7");

	MultiSequence InputSeqs;
	LoadInput(InputSeqs);

	const string &OutputFileName = opt(output);
	bool Nucleo = false;
	if (opt(nt))
		Nucleo = true;
	else if (opt(amino))
		Nucleo = false;
	else
		Nucleo = InputSeqs.GuessIsNucleo();

	TREEPERM TP = TP_None;
	if (optset_perm)
		TP = StrToTREEPERM(opt(perm));

	SetAlpha(Nucleo ? ALPHA_Nucleo : ALPHA_Amino);
	InitProbcons();

	Super4 S4;
	S4.Run(InputSeqs, TP);
	S4.m_FinalMSA.WriteMFA(OutputFileName);
	}
