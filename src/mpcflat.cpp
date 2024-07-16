#include "muscle.h"
#include "mpcflat.h"
#include "tree.h"
#include "locallock.h"

void MPCFlat::Clear()
	{
	m_MyInputSeqs = 0;
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
	m_Weights.clear();
	FreeSparsePosts();
	FreeProgMSAs();
	}

uint MPCFlat::GetMyInputSeqIndex(const string &Label) const
	{
	unordered_map<string, uint>::const_iterator iter =
	  m_LabelToIndex.find(Label);
	asserta(iter != m_LabelToIndex.end());
	uint Index = iter->second;
	return Index;
	}

const char *MPCFlat::GetLabel(uint SeqIndex) const
	{
	const char *Label = m_MyInputSeqs->GetLabel(SeqIndex);
	return Label;
	}

uint MPCFlat::GetSeqLength(uint SeqIndex) const
	{
	uint L = m_MyInputSeqs->GetSeqLength(SeqIndex);
	return L;
	}

const Sequence *MPCFlat::GetSequence(uint SeqIndex) const
	{
	const Sequence *s = m_MyInputSeqs->GetSequence(SeqIndex);
	return s;
	}

const byte *MPCFlat::GetBytePtr(uint SeqIndex) const
	{
	const byte *Ptr = m_MyInputSeqs->GetBytePtr(SeqIndex);
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

uint MPCFlat::GetSeqCount() const
	{
	asserta(m_MyInputSeqs != 0);
	uint SeqCount = m_MyInputSeqs->GetSeqCount();
	return SeqCount;
	}

void MPCFlat::InitSeqs(MultiSequence *InputSeqs)
	{
	m_MyInputSeqs = InputSeqs;
	const uint SeqCount = GetSeqCount();

	m_Labels.clear();
	m_LabelToIndex.clear();
	for (uint i = 0; i < SeqCount; ++i)
		{
		const Sequence *Seq = InputSeqs->GetSequence(i);

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
	if (optset_guidetreeout)
		{
		const string &FN = opt(guidetreeout);
		Progress("Saving guide tree [%s] to %s\n",
		  TREEPERMToStr(m_TreePerm), FN.c_str());
		m_GuideTree.ToFile(FN);
		Progress("Quitting.\n");
		exit(0);
		}
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

void MPCFlat::Run(MultiSequence *OriginalInputSeqs)
	{
	assert(OriginalInputSeqs != 0);
	Clear();

	m_D.Run(*OriginalInputSeqs);
	m_D.Validate();
	
	MultiSequence *UniqueSeqs = new MultiSequence;
	m_D.GetUniqueSeqs(*UniqueSeqs);

	const uint SeqCount = UniqueSeqs->GetSeqCount();
	if (SeqCount == 1)
		{
		m_MSA = new MultiSequence;
		m_MSA->Copy(*OriginalInputSeqs);
		return;
		}
	unordered_map<string, vector<string> > RepSeqLabelToDupeLabels;
	m_D.GetRepLabelToDupeLabels(RepSeqLabelToDupeLabels);
	const uint DupeCount = SIZE(RepSeqLabelToDupeLabels);

	const uint PairCount = (SeqCount*(SeqCount-1))/2;
	AllocPairCount(PairCount);

	InitSeqs(UniqueSeqs);
	InitPairs();
	InitDistMx();
	CalcPosteriors();
	CalcGuideTree();

	m_CW.Run(*m_MyInputSeqs, m_GuideTree, m_Weights);
	asserta(SIZE(m_Weights) == SeqCount);
	double Sum = 0;
	for (uint i = 0; i < SeqCount; ++i)
		{
		float w = m_Weights[i]*SeqCount;
		m_Weights[i] = w;
		Sum += w;
		m_Weights[i] = 1.0f;//@@@@@@@@@ TODO
		}
	asserta(feq(Sum, SeqCount));

	Consistency();
	CalcJoinOrder();
	ProgressiveAlign();
	Refine();
	asserta(m_MSA != 0);
	SortMSA();

	if (DupeCount > 0)
		InsertDupes(RepSeqLabelToDupeLabels);
	}

void MPCFlat::SortMSA()
	{
	if (opt(input_order))
		SortMSA_ByInputOrder();
	else
		{
		opt(tree_order);
		SortMSA_ByGuideTree();
		}
	}

void MPCFlat::GetLabelToMSASeqIndex(unordered_map<string, uint> &LabelToMSASeqIndex) const
	{
	LabelToMSASeqIndex.clear();
	const uint SeqCount = GetSeqCount();
	asserta(m_MSA->GetSeqCount() == SeqCount);
	for (uint i = 0; i < SeqCount; ++i)
		{
		const string &Label = m_MSA->GetLabelStr(i);
		if (LabelToMSASeqIndex.find(Label) != LabelToMSASeqIndex.end())
			Die("Duplicate label >%s", Label.c_str());
		LabelToMSASeqIndex[Label] = i;
		}
	}

void MPCFlat::SortMSA_ByGuideTree()
	{
	const uint SeqCount = GetSeqCount();
	unordered_map<string, uint> LabelToMSASeqIndex;
	GetLabelToMSASeqIndex(LabelToMSASeqIndex);

	const uint LeafCount = m_GuideTree.GetLeafCount();
	asserta(LeafCount == SeqCount);
	const Tree &T = m_GuideTree;
	uint Node = T.FirstDepthFirstNode();
	vector<const Sequence *> SortedSeqs;
	for (;;)
		{
		if (T.IsLeaf(Node))
			{
			string Label;
			T.GetLabel(Node, Label);
			unordered_map<string, uint>::const_iterator p =
			  LabelToMSASeqIndex.find(Label);
			if (p == LabelToMSASeqIndex.end())
				Die("MPCFlat::SortMSA_ByGuideTree not found >%s", Label.c_str());
			uint MSASeqIndex = p->second;
			const Sequence *Seq = m_MSA->GetSequence(MSASeqIndex);
			SortedSeqs.push_back(Seq);
			}

		Node = T.NextDepthFirstNode(Node);
		if (Node == NULL_NEIGHBOR)
			break;
		}

	for (uint i = 0; i < SeqCount; ++i)
		m_MSA->m_Seqs[i] = SortedSeqs[i];
	}

void MPCFlat::SortMSA_ByInputOrder()
	{
	const uint SeqCount = GetSeqCount();
	unordered_map<string, uint> LabelToMSASeqIndex;
	GetLabelToMSASeqIndex(LabelToMSASeqIndex);

	vector<const Sequence *> SortedSeqs;
	for (uint i = 0; i < SeqCount; ++i)
		{
		const string &Label = m_MyInputSeqs->GetLabelStr(i);
		unordered_map<string, uint>::const_iterator p =
		  LabelToMSASeqIndex.find(Label);
		if (p == LabelToMSASeqIndex.end())
			Die("MPCFlat::SortMSA_ByInputOrder(), missing >%s", Label.c_str());
		uint MSASeqIndex = p->second;
		const Sequence *Seq = m_MSA->GetSequence(MSASeqIndex);
		SortedSeqs.push_back(Seq);
		}
	for (uint i = 0; i < SeqCount; ++i)
		m_MSA->m_Seqs[i] = SortedSeqs[i];
	}

void MPCFlat::InsertDupes(
  const unordered_map<string, vector<string> > &RepSeqLabelToDupeLabels)
	{
	MultiSequence *UpdatedMSA = new MultiSequence;
	const uint SeqCount = m_MSA->GetSeqCount();
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const Sequence *OldSeq = m_MSA->GetSequence(SeqIndex);
		Sequence *NewSeq = Sequence::NewSequence();
		const string &Label = OldSeq->m_Label;
		NewSeq->m_Label = Label;
		NewSeq->m_CharVec = OldSeq->m_CharVec;

		UpdatedMSA->AddSequence(NewSeq, false);
		unordered_map<string, vector<string> >::const_iterator iter =
		  RepSeqLabelToDupeLabels.find(Label);
		if (iter != RepSeqLabelToDupeLabels.end())
			{
			const vector<string> &DupeLabels = iter->second;
			for (uint i = 0; i < SIZE(DupeLabels); ++i)
				{
				Sequence *DupeSeq = Sequence::_NewSequence();
				DupeSeq->m_Label = DupeLabels[i];
				DupeSeq->m_CharVec = OldSeq->m_CharVec;
				UpdatedMSA->AddSequence(DupeSeq, true);
				}
			}
		}
	delete m_MSA;
	m_MSA = UpdatedMSA;
	}
