#include "muscle.h"
#include "mpcflat.h"
#include "tree.h"
#include "locallock.h"

#define DOTIMING	0

#if DOTIMING
#include "timing.h"

static TICKS g_tFwd;
static TICKS g_tBwd;
static TICKS g_tPost;
static TICKS g_tSparse;
static TICKS g_tAln;
static TICKS g_tDelete;

#endif

void MPCFlat::Clear()
	{
	m_InputSeqs = 0;
	if (m_MSA != 0)
		delete m_MSA;
	m_MSA = 0;
	m_Labels.clear();
	m_LabelToIndex.clear();
	m_Upgma5.Clear();
	m_GuideTree.Clear();
	m_DistMx.clear();
	m_Pairs.clear();
	m_PairToIndex.clear();
	m_JoinIndexes1.clear();
	m_JoinIndexes2.clear();
	FreeSparsePosts();
	FreeProgMSAs();
	}


const char *MPCFlat::GetLabel(uint SeqIndex) const
	{
	const char *Label = m_InputSeqs->GetLabel(SeqIndex);
	return Label;
	}

uint MPCFlat::GetSeqLength(uint SeqIndex) const
	{
	uint L = m_InputSeqs->GetSeqLength(SeqIndex);
	return L;
	}

const Sequence *MPCFlat::GetSequence(uint SeqIndex) const
	{
	const Sequence *s = m_InputSeqs->GetSequence(SeqIndex);
	return s;
	}

const byte *MPCFlat::GetBytePtr(uint SeqIndex) const
	{
	const byte *Ptr = m_InputSeqs->GetBytePtr(SeqIndex);
	return Ptr;
	}

uint MPCFlat::GetPairIndex(uint SMI1, uint SMI2) const
	{
	//asserta(SMI1 > SMI2);
	//uint PairIndex = SMI2 + (SMI1*(SMI1 - 1))/2;
	asserta(SMI1 < SMI2);
	const pair<uint, uint> Pair(SMI1, SMI2);
	map<pair<uint, uint>, uint>::const_iterator p = m_PairToIndex.find(Pair);
	asserta(p != m_PairToIndex.end());
	uint PairIndex = p->second;
	return PairIndex;
	}

const pair<uint, uint> &MPCFlat::GetPair(uint PairIndex) const
	{
	assert(PairIndex < SIZE(m_Pairs));
	return m_Pairs[PairIndex];
	}

void MPCFlat::AllocPairCount(uint PairCount)
	{
	asserta(PairCount > 0);

	if (PairCount < SIZE(*m_ptrSparsePosts))
		return;

	m_SparsePosts1.resize(PairCount);
	m_SparsePosts2.resize(PairCount);
	}

MySparseMx &MPCFlat::GetSparsePost(uint PairIndex)
	{
	asserta(PairIndex < SIZE(*m_ptrSparsePosts));
	MySparseMx *Mx = (*m_ptrSparsePosts)[PairIndex];
	if (Mx == 0)
		{
		Mx = new MySparseMx;
		(*m_ptrSparsePosts)[PairIndex] = Mx;
		}
	return *Mx;
	}

MySparseMx &MPCFlat::GetUpdatedSparsePost(uint PairIndex)
	{
	asserta(PairIndex < SIZE(*m_ptrUpdatedSparsePosts));
	MySparseMx *Mx = (*m_ptrUpdatedSparsePosts)[PairIndex];
	if (Mx == 0)
		{
		Mx = new MySparseMx;
		(*m_ptrUpdatedSparsePosts)[PairIndex] = Mx;
		}
	return *Mx;
	}

uint MPCFlat::GetL(uint SeqIndex) const
	{
	return m_InputSeqs->GetSeqLength(SeqIndex);
	}

uint MPCFlat::GetSeqCount() const
	{
	asserta(m_InputSeqs != 0);
	uint SeqCount = m_InputSeqs->GetSeqCount();
	return SeqCount;
	}

void MPCFlat::InitSeqs(MultiSequence *InputSeqs)
	{
	m_InputSeqs = InputSeqs;
	const uint SeqCount = GetSeqCount();

	m_Labels.clear();
	m_LabelToIndex.clear();
	for (uint i = 0; i < SeqCount; ++i)
		{
		const Sequence *Seq = InputSeqs->GetSequence(i);
		Sequence *HackSeq = (Sequence *) Seq;
		HackSeq->m_SMI = i;

		const string &Label = Seq->GetLabel();
		m_Labels.push_back(Label);

		if (m_LabelToIndex.find(Label) != m_LabelToIndex.end())
			Die("Duplicate label >%s", Label.c_str());
		m_LabelToIndex[Label] = i;
		}
	}

void MPCFlat::InitPairs()
	{
	const uint SeqCount = GetSeqCount();
	m_Pairs.clear();
	m_PairToIndex.clear();
	uint PairIndex = 0;
	for (uint SeqIndex1 = 0; SeqIndex1 < SeqCount; ++SeqIndex1)
		for (uint SeqIndex2 = SeqIndex1 + 1; SeqIndex2 < SeqCount; ++SeqIndex2)
			{
			const pair<uint, uint> Pair(SeqIndex1, SeqIndex2);
			m_Pairs.push_back(Pair);
			assert(m_PairToIndex.find(Pair) == m_PairToIndex.end());
			m_PairToIndex[Pair] = PairIndex;
			uint PairIndex2 = GetPairIndex(SeqIndex1, SeqIndex2);
			asserta(PairIndex2 == PairIndex);
			++PairIndex;
			}
	uint PairCount = (SeqCount * (SeqCount - 1)) / 2;
	uint PairCount2 = SIZE(m_Pairs);
	asserta(PairCount == PairCount2);
	}

void MPCFlat::InitDistMx()
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

void MPCFlat::Consistency()
	{
	const uint SeqCount = GetSeqCount();
	if (SeqCount < 3)
		return;

	for (uint Iter = 0; Iter < m_ConsistencyIterCount; ++Iter)
		ConsIter(Iter);
	}

void MPCFlat::CalcGuideTree()
	{
	if (opt(randomchaintree))
		{
		CalcGuideTree_RandomChain();
		return;
		}

	m_Upgma5.Init(m_Labels, m_DistMx);
	m_Upgma5.FixEADistMx();

	m_Upgma5.Run(LINKAGE_Biased, m_GuideTree);
	PermTree(m_GuideTree, m_TreePerm);
	}

void MPCFlat::CalcJoinOrder()
	{
	GetGuideTreeJoinOrder(m_GuideTree, m_LabelToIndex,
	  m_JoinIndexes1, m_JoinIndexes2);
	ValidateJoinOrder(m_JoinIndexes1, m_JoinIndexes2);
	}

void MPCFlat::CalcPosteriors()
	{
#if 0//TRACE
	{
	Log("MPCFlat::CalcPosteriors SIZE(m_SparsePosts1)=%u SIZE(m_SparsePosts2)=%u\n",
	  SIZE(m_SparsePosts1), SIZE(m_SparsePosts2));
	for (uint i = 0; i < SIZE(m_SparsePosts1); ++i)
		{
		const MySparseMx *M = m_SparsePosts1[i];
		Log("m_SparsePosts1[%u]=[%p]", i, M);
		if (M != 0)
			Log(" maxlx %u\n", M->m_MaxLX);
		Log("\n");
		}
	for (uint i = 0; i < SIZE(m_SparsePosts2); ++i)
		{
		const MySparseMx *M = m_SparsePosts2[i];
		Log("m_SparsePosts2[%u]=[%p]", i, M);
		if (M != 0)
			Log(" maxlx %u\n", M->m_MaxLX);
		Log("\n");
		}
	}
#endif

	uint PairCount = SIZE(m_Pairs);
	asserta(PairCount > 0);
	unsigned ThreadCount = GetRequestedThreadCount();
	uint PairCounter = 0;
#pragma omp parallel for num_threads(ThreadCount)
	for (int PairIndex = 0; PairIndex < (int) PairCount; ++PairIndex)
		{
		Lock();
		ProgressStep(PairCounter++, PairCount, "Calc posteriors");
		Unlock();

		CalcPosterior(PairIndex);
		}
	}

void MPCFlat::Refine()
	{
	const uint SeqCount = GetSeqCount();
	if (SeqCount < 3)
		return;
	for (uint Iter = 0; Iter < m_RefineIterCount; ++Iter)
		{
		ProgressStep(Iter, m_RefineIterCount, "Refining");
		RefineIter();
		}
	}

void MPCFlat::Run_Super4(MultiSequence *ConsensusSeqs)
	{
	assert(ConsensusSeqs != 0);
	Clear();

	const uint SeqCount = ConsensusSeqs->GetSeqCount();
	asserta(SeqCount > 1);

	uint PairCount = (SeqCount*(SeqCount-1))/2;
	AllocPairCount(PairCount);

	InitSeqs(ConsensusSeqs);
	InitPairs();
	InitDistMx();
	CalcPosteriors();
	Consistency();
	CalcGuideTree();
	}

void MPCFlat::Run(MultiSequence *InputSeqs)
	{
	assert(InputSeqs != 0);
	Clear();

	const uint SeqCount = InputSeqs->GetSeqCount();
	if (SeqCount == 1)
		{
		m_MSA = InputSeqs;
		return;
		}

	uint PairCount = (SeqCount*(SeqCount-1))/2;
	AllocPairCount(PairCount);

	InitSeqs(InputSeqs);
	InitPairs();
	InitDistMx();
	CalcPosteriors();
	Consistency();
	CalcGuideTree();
	CalcJoinOrder();
	ProgressiveAlign();
	Refine();
	asserta(m_MSA != 0);
	}

MultiSequence *RunMPCFlat(MultiSequence *InputSeqs)
	{
	MPCFlat M;
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

	asserta(M.m_MSA != 0);
	return M.m_MSA;
	}

#if DOTIMING
void LogTiming()
	{
	double Sum = g_tFwd + g_tBwd + g_tPost + g_tAln + g_tDelete;
	ProgressLog("%10.3g  %4.1f%%  Fwd\n", double(g_tFwd), GetPct(g_tFwd, Sum));
	ProgressLog("%10.3g  %4.1f%%  Bwd\n", double(g_tBwd), GetPct(g_tBwd, Sum));
	ProgressLog("%10.3g  %4.1f%%  Post\n", double(g_tPost), GetPct(g_tPost, Sum));
	ProgressLog("%10.3g  %4.1f%%  Aln\n", double(g_tAln), GetPct(g_tAln, Sum));
	ProgressLog("%10.3g  %4.1f%%  Delete\n", double(g_tDelete), GetPct(g_tDelete, Sum));
	}
#endif
