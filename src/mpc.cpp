#include "muscle.h"
#include "tree.h"
#include "mpc.h"

uint MPC::GetSeqCount() const
	{
	asserta(m_InputSeqs != 0);
	uint SeqCount = m_InputSeqs->GetSeqCount();
	return SeqCount;
	}

void MPC::InitSeqs(MultiSequence *InputSeqs)
	{
	m_InputSeqs = InputSeqs;
	const uint SeqCount = GetSeqCount();

	m_Labels.clear();
	m_LabelToIndex.clear();
	for (uint i = 0; i < SeqCount; ++i)
		{
		Sequence *Seq = InputSeqs->GetSequence(i);
		Seq->SetSMI(i);

		const string &Label = Seq->GetLabel();
		m_Labels.push_back(Label);

		if (m_LabelToIndex.find(Label) != m_LabelToIndex.end())
			Die("Duplicate label >%s", Label.c_str());
		m_LabelToIndex[Label] = i;
		}
	}

void MPC::InitPairs()
	{
	const uint SeqCount = GetSeqCount();
	m_Pairs.clear();
	for (uint SeqIndex1 = 0; SeqIndex1 < SeqCount; ++SeqIndex1)
		for (uint SeqIndex2 = SeqIndex1 + 1; SeqIndex2 < SeqCount; ++SeqIndex2)
			m_Pairs.push_back(pair<uint, uint>(SeqIndex1, SeqIndex2));
	uint PairCount = (SeqCount * (SeqCount - 1)) / 2;
	uint PairCount2 = SIZE(m_Pairs);
	asserta(PairCount == PairCount2);
	}

void MPC::InitDistMx()
	{
	const uint SeqCount = GetSeqCount();
	m_DistMx.clear();
	m_DistMx.resize(SeqCount);
	for (uint i = 0; i < SeqCount; ++i)
		{
		m_DistMx[i].resize(SeqCount, FLT_MAX);
		m_DistMx[i][i] = 0;
		}
	}

void MPC::CalcPosterior(uint SeqIndex1, uint SeqIndex2)
	{
	Sequence* seq1 = m_InputSeqs->GetSequence(SeqIndex1);
	Sequence* seq2 = m_InputSeqs->GetSequence(SeqIndex2);

	uint L1 = seq1->GetLength();
	uint L2 = seq2->GetLength();

	vector<float>* forward = PairHMM::ComputeForwardMatrix(seq1, seq2);
	asserta(forward != 0);

	vector<float>* backward = PairHMM::ComputeBackwardMatrix(seq1, seq2);
	asserta(backward != 0);

	vector<float>* posterior =
		PairHMM::ComputePosteriorMatrix(seq1, seq2, *forward, *backward);
	asserta(posterior != 0);

	m_SparseMatrices[SeqIndex1][SeqIndex2] = new SparseMatrix(L1, L2, *posterior);
	m_SparseMatrices[SeqIndex2][SeqIndex1] = NULL;

	pair<vector<char>*, float> alignment =
		PairHMM::ComputeAlignment(L1, L2, *posterior);

// "expected accuracy" distance
	float EA = alignment.second / min(L1, L2);
	m_DistMx[SeqIndex1][SeqIndex2] = EA;
	m_DistMx[SeqIndex2][SeqIndex1] = EA;

	delete alignment.first;
	delete posterior;
	delete forward;
	delete backward;
	}

void MPC::Consistency()
	{
	const uint SeqCount = GetSeqCount();
	if (SeqCount < 3)
		return;

	for (uint Iter = 0; Iter < m_ConsistencyIterCount; ++Iter)
		DoRelax(Iter);
	}

void MPC::InitSparseMatrices()
	{
	const uint SeqCount = GetSeqCount();
	m_SparseMatrices.clear();
	m_SparseMatrices.resize(SeqCount);
	for (uint i = 0; i < SeqCount; ++i)
		m_SparseMatrices[i].resize(SeqCount, 0);
	}

void MPC::CalcGuideTree()
	{
	m_U5.Init(m_Labels, m_DistMx);
	m_U5.FixEADistMx();

	m_U5.Run(LINKAGE_Biased, m_GuideTree);
	PermTree(m_GuideTree, m_TreePerm);
	}

void MPC::CalcJoinOrder()
	{
	GetGuideTreeJoinOrder(m_GuideTree, m_LabelToIndex,
	  m_JoinIndexes1, m_JoinIndexes2);
	ValidateJoinOrder(m_JoinIndexes1, m_JoinIndexes2);
	}

void MPC::Refine()
	{
	const uint SeqCount = GetSeqCount();
	if (SeqCount < 3)
		return;

	for (uint RefineIter = 0; RefineIter < m_RefineIterCount; RefineIter++)
		{
		ProgressStep(RefineIter, m_RefineIterCount, "Refining");
		DoIterativeRefinement(m_SparseMatrices, m_MSA);
		}
	}

void MPC::CalcPosteriors()
	{
	const uint SeqCount = GetSeqCount();
	uint PairIndex = 0;
	uint PairCount = SIZE(m_Pairs);
	asserta(PairCount > 0);
	unsigned ThreadCount = GetRequestedThreadCount();
	omp_lock_t Lock;
	omp_init_lock(&Lock);
	uint PairCounter = 0;
#pragma omp parallel for num_threads(ThreadCount)
	for (int PairIndex = 0; PairIndex < (int) PairCount; ++PairIndex)
		{
		const pair<uint, uint>& Pair = m_Pairs[PairIndex];
		uint SeqIndex1 = Pair.first;
		uint SeqIndex2 = Pair.second;

		omp_set_lock(&Lock);
		ProgressStep(PairCounter++, PairCount, "Calc posteriors");
		omp_unset_lock(&Lock);

		CalcPosterior(SeqIndex1, SeqIndex2);
		}
	}

void MPC::DeleteSparseMatrices()
	{
	const uint SeqCount = GetSeqCount();
	for (uint SeqIndex1 = 0; SeqIndex1 < SeqCount - 1; SeqIndex1++)
		{
		for (uint SeqIndex2 = SeqIndex1 + 1; SeqIndex2 < SeqCount; SeqIndex2++)
			{
			delete m_SparseMatrices[SeqIndex1][SeqIndex2];
			delete m_SparseMatrices[SeqIndex2][SeqIndex1];
			}
		}
	m_SparseMatrices.clear();
	}

void MPC::Run(MultiSequence *InputSeqs)
	{
	assert(InputSeqs != 0);

	const uint SeqCount = InputSeqs->GetSeqCount();
	if (SeqCount == 1)
		{
		m_MSA = InputSeqs;
		return;
		}

	InitSeqs(InputSeqs);
	InitPairs();
	InitDistMx();
	InitSparseMatrices();
	CalcPosteriors();
	Consistency();
	CalcGuideTree();
	CalcJoinOrder();
	ProgressiveAlign();
	Refine();
	CalcAccAlnRows();
	DeleteSparseMatrices();
	asserta(m_MSA != 0);
	}

MultiSequence *RunMPC(MultiSequence *InputSeqs)
	{
	MPC M;
	if (optset_consiters)
		M.m_ConsistencyIterCount = opt(consiters);
	if (optset_refineiters)
		M.m_RefineIterCount = opt(refineiters);

	TREEPERM TP = TP_None;
	if (optset_perm)
		TP = StrToTREEPERM(opt(perm));
	if (TP == TP_All)
		Die("-perm all not supported, please specify none, abc, acb or bca");
	M.m_TreePerm = TP;
	M.Run(InputSeqs);
	M.WriteAccAln(opt(accalnout));

	asserta(M.m_MSA != 0);
	return M.m_MSA;
	}
