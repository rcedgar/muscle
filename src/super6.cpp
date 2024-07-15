#include "muscle.h"
#include "super6.h"
#include "omplock.h"

double GetProtDistMFAPair(const MultiSequence &MFA1,
  const MultiSequence &MFA2, uint TargetPairCount);

void Super6::SetOpts()
	{
	m_MaxPDPass1 = (float) optd(super6_maxpd1, DEFAULT_MAX_PD_SUPER6_PASS1);
	}

void Super6::Run(MultiSequence &InputSeqs)
	{
	m_InputSeqs = &InputSeqs;
	const uint SeqCount = InputSeqs.GetSeqCount();
	vector<uint> AllSeqIndexes;
	for (uint i = 0; i < SeqCount; ++i)
		AllSeqIndexes.push_back(i);

	m_UCPD.Run(InputSeqs, AllSeqIndexes, m_MaxPDPass1);
	m_UCPD.GetClusterMFAs(m_ClusterMFAs);
	AssertSameSeqsVec(*m_InputSeqs, m_ClusterMFAs);
	PrepareClusters();
	CalcClusterDistMx();
	MakeGuideTree();
	AlignClusters();
	InitPP();
	m_PP.RunGuideTree(m_GuideTree);
	}

void Super6::AlignClusters()
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

void Super6::SplitBigMFA_Random(MultiSequence &InputMFA, uint MaxSize,
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

void Super6::PrepareClusters()
	{
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
			SplitBigMFA_Random(*MFA, m_MaxClusterSize, *SplitMFAs);

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

void Super6::CalcClusterDistMx()
	{
	const uint ClusterCount = SIZE(m_ClusterMFAs);
	m_ClusterDistMx.clear();
	m_ClusterDistMx.resize(ClusterCount);
	for (uint i = 0; i < ClusterCount; ++i)
		{
		m_ClusterDistMx[i].resize(ClusterCount, FLT_MAX);
		m_ClusterDistMx[i][i] = 0;
		}

	const uint PairCount = (ClusterCount*(ClusterCount - 1))/2;
	unsigned ThreadCount = GetRequestedThreadCount();
	uint ClusterIndex1 = 1;
	uint ClusterIndex2 = 0;
	uint PairCounter = 0;
#pragma omp parallel for num_threads(ThreadCount)
	for (int PairIndex = 0; PairIndex < (int) PairCount; ++PairIndex)
		{
		LOCK();
		ProgressStep(PairCounter++, PairCount, "Cluster dist mx");
		UNLOCK();

		const MultiSequence &MFA1 = *m_ClusterMFAs[ClusterIndex1];
		const MultiSequence &MFA2 = *m_ClusterMFAs[ClusterIndex2];
		float d = (float) GetProtDistMFAPair(MFA1, MFA2, 8);
		m_ClusterDistMx[ClusterIndex1][ClusterIndex2] = d;
		m_ClusterDistMx[ClusterIndex2][ClusterIndex1] = d;

		LOCK();
		++ClusterIndex2;
		if (ClusterIndex2 == ClusterIndex1)
			{
			++ClusterIndex1;
			ClusterIndex2 = 0;
			}
		UNLOCK();
		}
	}

void Super6::MakeGuideTree()
	{
	UPGMA5 U;
	U.Init(m_ClusterLabels, m_ClusterDistMx);
	U.Run(LINKAGE_Biased, m_GuideTree);
	m_GuideTree.LogMe();
	m_GuideTree.ToFile("guide.newick");
	}

void Super6::InitPP()
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

void cmd_super6()
	{
	MultiSequence InputSeqs;
	LoadInput(InputSeqs);

	//LoadGlobalInputMS(opt(super6));

	string &OutputPattern = opt(output);
	if (OutputPattern.empty())
		Die("Must set -output");

	//MultiSequence &InputSeqs = GetGlobalInputMS();
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
	SetAlphab(Nucleo);

	Super6 S6;
	S6.SetOpts();
	S6.Run(InputSeqs);
	const MultiSequence &FinalMSA = S6.m_PP.GetFinalMSA();
	FinalMSA.WriteMFA(opt(output));
	}
