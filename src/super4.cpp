#include "muscle.h"
#include "super4.h"
#include "sequence.h"
#include "multisequence.h"
#include "usorter.h"
#include "upgma5.h"
#include "pprog.h"
#include "treeperm.h"

TREEPERM StrToTREEPERM(const string &s)
	{
	if (s == "none") return TP_None;
	if (s == "abc") return TP_ABC;
	if (s == "acb") return TP_ACB;
	if (s == "bca") return TP_BCA;
	if (s == "all") return TP_All;
	Die("Invalid perm '%s'", s.c_str());
	return TP_None;
	}

const char *TREEPERMToStr(TREEPERM TP)
	{
	switch (TP)
		{
	case TP_None:	return "none";
	case TP_ABC:	return "abc";
	case TP_ACB:	return "acb";
	case TP_BCA:	return "bca";
	case TP_All:	return "all";
	default:
		break;
		}
	asserta(false);
	return "?";
	}

void Super4::MakeGuideTree()
	{
	UPGMA5 U;
	U.Init(m_ClusterLabels, m_DistMx);
	U.FixEADistMx();
	U.Run(LINKAGE_Biased, m_GuideTree);
	PermTree(m_GuideTree, m_TreePerm);
	}

void Super4::CalcDistMx()
	{
	FILE *f = 0;
	if (m_SaveDir != "")
		{
		string FileName;
		Ps(FileName, "%sdistmx.tsv", m_SaveDir.c_str());
		f = CreateStdioFile(FileName);
		}

	CalcEADistMx(f, &m_ConsensusSeqs, m_DistMx);
	CloseStdioFile(f);
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
			Sequence *seq = InputMFA.GetSequence(OutputSeqCount + i);
			SplitMFA->AddSequence(seq);
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

void Super4::ClusterInput(uint MaxClusterSize)
	{
	const float MIN_EA_PASS1 = 0.7f;
	const float MIN_EA_PASS2 = 0.9f;
	const uint MAX_CLUSTER_SIZE = 500;

	const uint InputSeqCount = m_InputSeqs->GetSeqCount();

	m_EC.Run(*m_InputSeqs, MIN_EA_PASS1);
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
		if (SeqCount > MAX_CLUSTER_SIZE)
			{
			vector<MultiSequence *> *SplitMFAs = new vector<MultiSequence *>;
			SplitBigMFA(*MFA, MAX_CLUSTER_SIZE, MIN_EA_PASS2, *SplitMFAs);

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
		asserta(SeqCount <= MaxClusterSize);

		string ClusterLabel;
		Ps(ClusterLabel, "Cluster%u", ClusterIndex);
		m_ClusterLabels.push_back(ClusterLabel);
		}
	}

void Super4::WriteMSAs() const
	{
	if (m_SaveDir == "")
		return;

	uint CentroidCount = SIZE(m_ClusterMSAs);
	for (uint CentroidIndex = 0; CentroidIndex < CentroidCount; ++CentroidIndex)
		{
		ProgressStep(CentroidIndex, CentroidCount, "Write MSAs %s",
		  m_SaveDir.c_str());
		const MultiSequence *MSA = m_ClusterMSAs[CentroidIndex];
		string FileName;
		Ps(FileName, "%sclust%u.afa", m_SaveDir.c_str(), CentroidIndex);
		MSA->WriteMFA(FileName);
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

		AssertSameLabels(*m_InputSeqs);
		MultiSequence *ClusterMSA = RunMPC(ClusterMFA);
		AssertSameLabels(*m_InputSeqs);
		AssertSameLabels(*ClusterMSA);
		m_ClusterMSAs.push_back(ClusterMSA);
		}
	}

void Super4::GetConsensusSeqs()
	{
	uint ClusterCount = SIZE(m_ClusterMSAs);
	for (uint ClusterIndex = 0; ClusterIndex < ClusterCount; ++ClusterIndex)
		{
		ProgressStep(ClusterIndex, ClusterCount, "Consensus sequences");
		MultiSequence *ClusterMSA = m_ClusterMSAs[ClusterIndex];

		string Label;
		Ps(Label, "Cluster%u", ClusterIndex);

		string ConsSeq;
		GetConsensusSequence(*ClusterMSA, ConsSeq);

		Sequence *seq = new Sequence;
		seq->FromString(Label, ConsSeq);

		m_ConsensusSeqs.AddSequence(seq);

		asserta(ClusterIndex < SIZE(m_ClusterLabels));
		const string &ClusterLabel = m_ClusterLabels[ClusterIndex];
		}
	}

void Super4::WriteConsensusSeqs()
	{
	if (m_SaveDir == "")
		return;

	string FileName;
	Ps(FileName, "%scons.fa", m_SaveDir.c_str());
	m_ConsensusSeqs.WriteMFA(FileName);
	}

void Super4::Final()
	{
	vector<const MultiSequence *> MSAs;
	const uint n = SIZE(m_ClusterMSAs);
	for (uint i = 0; i < n; ++i)
		{
		const MultiSequence *MSA = m_ClusterMSAs[i];
		MSAs.push_back(MSA);
		}

	m_PP.m_TargetPairCount = m_PairCount;
	m_PP.SetMSAs(MSAs, m_ClusterLabels);
	m_PP.RunGuideTree(m_GuideTree);
	const MultiSequence &FinalMSA = m_PP.GetFinalMSA();
	m_FinalMSA.Copy(FinalMSA);

	//PermuteTree(m_GuideTree, m_GuideTreeABC, m_GuideTreeACB,
	//  m_GuideTreeBCA, m_LabelsA, m_LabelsB, m_LabelsC);

	//if (m_SaveDir != "")
	//	{
	//	m_GuideTreeABC.ToFile(m_SaveDir + "ABC.newick");
	//	m_GuideTreeACB.ToFile(m_SaveDir + "ACB.newick");
	//	m_GuideTreeBCA.ToFile(m_SaveDir + "BCA.newick");

	//	StringsToFile(m_SaveDir + "labelsA.txt", m_LabelsA);
	//	StringsToFile(m_SaveDir + "labelsB.txt", m_LabelsB);
	//	StringsToFile(m_SaveDir + "labelsC.txt", m_LabelsC);
	//	}

	//m_PP_ABC.m_TargetPairCount = m_PairCount;
	//m_PP_ABC.SetMSAs(MSAs, m_ClusterLabels);
	//m_PP_ABC.RunGuideTree(m_GuideTreeABC);
	//m_FinalMSA_ABC.Copy(m_PP_ABC.GetFinalMSA());

	//m_PP_ACB.m_TargetPairCount = m_PairCount;
	//m_PP_ACB.SetMSAs(MSAs, m_ClusterLabels);
	//m_PP_ACB.RunGuideTree(m_GuideTreeACB);
	//m_FinalMSA_ACB.Copy(m_PP_ACB.GetFinalMSA());

	//m_PP_BCA.m_TargetPairCount = m_PairCount;
	//m_PP_BCA.SetMSAs(MSAs, m_ClusterLabels);
	//m_PP_BCA.RunGuideTree(m_GuideTreeBCA);
	//m_FinalMSA_BCA.Copy(m_PP_BCA.GetFinalMSA());

	//if (m_SaveDir != "")
	//	{
	//	m_FinalMSA_ABC.WriteMFA(m_SaveDir + "ABC.afa");
	//	m_FinalMSA_ACB.WriteMFA(m_SaveDir + "ACB.afa");
	//	m_FinalMSA_BCA.WriteMFA(m_SaveDir + "BCA.afa");
	//	}
	}

void RunSuper4(MultiSequence &InputSeqs, MultiSequence &MSA, TREEPERM TreePerm)
	{
	const uint PairCount = optd(paircount, DEFAULT_TARGET_PAIR_COUNT);
	const uint MaxCoarse = optd(maxcoarse, 500);
	string SaveDir;
	if (optset_savedir)
		{
		SaveDir = opt(savedir);
		Dirize(SaveDir);
		}
	asserta(!optset_minea);
	AssertSameLabels(InputSeqs);

	SetAlpha(ALPHA_Amino);

	Super4 S4;
	S4.m_SaveDir = SaveDir;
	InitProbcons();
	S4.m_TreePerm = TreePerm;
	S4.m_InputSeqs = &InputSeqs;
	S4.m_fTSV = 0;
	if (S4.m_SaveDir != "")
		S4.m_fTSV = CreateStdioFile(S4.m_SaveDir + "super4.tsv");

	S4.m_PairCount = PairCount;

	S4.ClusterInput(MaxCoarse);
	S4.AlignClusters();
	S4.WriteMSAs();
	S4.GetConsensusSeqs();
	S4.WriteConsensusSeqs();
	S4.CalcDistMx();
	S4.MakeGuideTree();
	S4.Final();
	MSA.Copy(S4.m_FinalMSA);

	AssertSameSeqs(InputSeqs, MSA);

	CloseStdioFile(S4.m_fTSV);
	}

void cmd_super4()
	{
	const string &InputFileName = opt(super4);
	const string &OutputFileName = opt(output);

	MultiSequence &InputSeqs = LoadGlobalInputMS(InputFileName);

	TREEPERM TP = TP_None;
	if (optset_perm)
		TP = StrToTREEPERM(opt(perm));
	if (TP == TP_All)
		Die("-perm all not supported, please specify none, abc, acb or bca");

	MultiSequence MSA;
	RunSuper4(InputSeqs, MSA, TP);
	MSA.WriteMFA(OutputFileName);
	}
