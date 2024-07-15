#include "muscle.h"
#include "kmerdist33.h"
#include "kmerdist66.h"
#include "clustalweights.h"
#include "pprog3.h"
#include "muscle3.h"

void Muscle3::Run(const M3AlnParams &AP, const MultiSequence &InputSeqs)
	{
	m_AP = &AP;
	asserta(m_AP->m_Ready);

	m_InputSeqs = &InputSeqs;
	const uint SeqCount = InputSeqs.GetSeqCount();

	const string &KD = m_AP->m_KmerDist;
	if (KD == "66")
		m_K66.GetDistMx(InputSeqs, m_DistMx);
	else if (KD == "33")
		m_K33.GetDistMx(InputSeqs, m_DistMx);
	else
		Die("Muscle3::Run, m_AP->m_KmerDist=%s", KD.c_str());

	m_Labels.resize(0);
	for (uint i = 0; i < SeqCount; ++i)
		{
		string Label = (string) InputSeqs.GetLabel(i);
		m_Labels.push_back(Label);
		}
#if DOSC
	m_SC.Run(m_DistMx, m_Labels, m_AP->m_Linkage, false);
	m_SC.GetTree(m_GuideTree);
#else
	m_U5.Init(m_Labels, m_DistMx);
	m_U5.Run(AP.m_Linkage, m_GuideTree);
#endif
	m_CW.Run(InputSeqs, m_GuideTree, m_InputSeqWeights);
	m_PP3.m_AP = m_AP;
	m_PP3.Run(InputSeqs, m_InputSeqWeights, m_GuideTree);
	assert(m_PP3.m_MSA.IsAligned());

	MultiSequence MSA_InputOrder;
	const uint TreeIters = AP.m_TreeIters;
	for (uint TreeIter = 0; TreeIter < TreeIters; ++TreeIter)
		{
		MSA_InputOrder.Clear();
		MSA_InputOrder.m_Seqs.resize(SeqCount, 0);
		MSA_InputOrder.m_Owners.resize(SeqCount, false);
		for (uint k = 0; k < SeqCount; ++k)
			{
			const Sequence *Seq = m_PP3.m_MSA.GetSequence(k);
			uint SeqIndex = GetGSIByLabel(Seq->m_Label);
			asserta(SeqIndex < SeqCount);
			MSA_InputOrder.m_Seqs[SeqIndex] = Seq;
			}

		GetKimuraDistMx(MSA_InputOrder, m_DistMx);
		m_AP->PerturbDistMx(m_DistMx);

#if DOSC
		m_SC.Run(m_DistMx, m_Labels, m_AP->m_Linkage, false);
		m_SC.GetTree(m_GuideTree);
#else
		m_U5.Init(m_Labels, m_DistMx);
		m_U5.Run(AP.m_Linkage, m_GuideTree);
#endif

		m_CW.Run(InputSeqs, m_GuideTree, m_InputSeqWeights);
		m_PP3.Run(InputSeqs, m_InputSeqWeights, m_GuideTree);
		assert(m_PP3.m_MSA.IsAligned());
		}
	m_FinalMSA = &m_PP3.m_MSA;
	}

void Muscle3::WriteMSA(const string &FileName) const
	{
	if (FileName == "")
		return;
	asserta(m_FinalMSA != 0);
	m_FinalMSA->WriteMFA(FileName);
	}

void Muscle3::RunRO(const M3AlnParams &AP, const MultiSequence &InputSeqs)
	{
	asserta(m_AP != 0);
	asserta(m_AP->m_Ready);
	m_InputSeqs = &InputSeqs;
	const uint SeqCount = InputSeqs.GetSeqCount();
	vector<uint> Order;
	for (uint i = 0; i < SeqCount; ++i)
		Order.push_back(i);
	Shuffle(Order);

	vector<float> Weights;
	Weights.push_back(1.0f);

	uint SeqIndex0 = Order[0];
	const Sequence &Seq0 = *InputSeqs.GetSequence(SeqIndex0);
	MultiSequence *AccumulatedMSA = new MultiSequence;
	AccumulatedMSA->AddSequence(&Seq0, false);
	Profile3 *AccumulatedProf = new Profile3;
	AccumulatedProf->FromMSA(*AccumulatedMSA,
	  m_AP->m_SubstMx_Letter, m_AP->m_GapOpen, Weights);

	CacheMem3 CM;
	for (uint k = 1; k < SeqCount; ++k)
		{
	// Profile for k'th sequence (size = 1)
		uint SeqIndexk = Order[k];
		const Sequence &Seqk = *InputSeqs.GetSequence(SeqIndexk);
		MultiSequence *MSAk = new MultiSequence;
		MSAk->AddSequence(&Seqk, false);
		Profile3 *Profk = new Profile3;
		Profk->FromMSA(*MSAk, m_AP->m_SubstMx_Letter, m_AP->m_GapOpen, Weights);

	// Align k'th sequence to accumulated profile
		string Path;
		NWSmall3(CM, *Profk, *AccumulatedProf, Path);

		float RelativeWeightOneSeq = 1.0f/k;
		float RelativeWeightAccumulatedSeqs = float(k);
		float Sum = RelativeWeightOneSeq + RelativeWeightAccumulatedSeqs;

		float w1 = RelativeWeightOneSeq/Sum;
		float wn = RelativeWeightAccumulatedSeqs/Sum;

		Profile3 *CombinedProf = new Profile3;
		AlignTwoProfsGivenPath(*Profk, w1, *AccumulatedProf, wn,
		  m_AP->m_SubstMx_Letter, m_AP->m_GapOpen,
		  Path, *CombinedProf);

		delete AccumulatedProf;
		AccumulatedProf = CombinedProf;

		MultiSequence *CombinedMSA = new MultiSequence;
		AlignTwoMSAsGivenPath(*MSAk, *AccumulatedMSA,
		  Path, *CombinedMSA);
		delete AccumulatedMSA;
		delete MSAk;
		AccumulatedMSA = CombinedMSA;
		}
	m_FinalMSA = AccumulatedMSA;
	}
