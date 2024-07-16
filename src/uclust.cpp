#include "muscle.h"
#include "uclust.h"
#include "pwpath.h"

float UClust::AlignSeqPair(uint SeqIndex1, uint SeqIndex2, string &Path)
	{
	Path.clear();
	const Sequence *Seq1 = m_InputSeqs->GetSequence(SeqIndex1);
	const Sequence *Seq2 = m_InputSeqs->GetSequence(SeqIndex2);
	const string &Label1 = Seq1->m_Label;
	const string &Label2 = Seq2->m_Label;
	asserta(Seq1 != 0);
	asserta(Seq2 != 0);
	float EA = AlignPairFlat(Label1, Label2, Path);
	return EA;
	}

void UClust::AddSeqToIndex(uint SeqIndex)
	{
	const Sequence *Seq = m_InputSeqs->GetSequence(SeqIndex);
	const byte *ByteSeq = Seq->GetBytePtr();
	const uint L = Seq->GetLength();
	m_US.AddSeq(ByteSeq, L, SeqIndex);
	}

uint UClust::Search(uint SeqIndex, string &Path)
	{
	const Sequence *Seq = m_InputSeqs->GetSequence(SeqIndex);
	const byte *ByteSeq = Seq->GetBytePtr();
	const uint L = Seq->GetLength();

	vector<uint> TopSeqIndexes;
	vector<uint> TopWordCounts;
	m_US.SearchSeq(ByteSeq, L, TopSeqIndexes, TopWordCounts);
	uint TopCount = SIZE(TopSeqIndexes);
	asserta(SIZE(TopWordCounts) == TopCount);
	if (TopCount == 0)
		return UINT_MAX;
	if (TopCount > MAX_REJECTS)
		TopCount = MAX_REJECTS;
	uint ThreadCount = GetRequestedThreadCount();

	uint CentroidSeqIndex = UINT_MAX;
	for (int TopIndex = 0; TopIndex < (int) TopCount; ++TopIndex)
		{
		uint TopSeqIndex = TopSeqIndexes[TopIndex];
		float EA = AlignSeqPair(SeqIndex, TopSeqIndex, Path);
		if (EA >= m_MinEA)
			{
			CentroidSeqIndex = TopSeqIndex;
			break;
			}
		}
	return CentroidSeqIndex;
	}

void UClust::Run(MultiSequence &InputSeqs, float MinEA)
	{
	m_InputSeqs = &InputSeqs;
	m_MinEA = MinEA;
	m_US.Init();

	const uint InputSeqCount = m_InputSeqs->GetSeqCount();
	const uint GSICount = GetGlobalMSSeqCount();
	vector<uint> GSIToInputSeqIndex(GSICount, UINT_MAX);
	for (uint SeqIndex = 0; SeqIndex < InputSeqCount; ++SeqIndex)
		{
		const string Label = string(m_InputSeqs->GetLabel(SeqIndex));
		uint GSI = GetGSIByLabel(Label);
		asserta(GSI < GSICount);
		asserta(GSIToInputSeqIndex[GSI] == UINT_MAX);
		GSIToInputSeqIndex[GSI] = SeqIndex;
		}

	uint CentroidCount = 0;
	uint MemberCount = 0;

	m_SeqIndexToCentroidSeqIndex.clear();
	m_SeqIndexToPath.clear();

	m_SeqIndexToCentroidSeqIndex.resize(InputSeqCount, UINT_MAX);
	m_SeqIndexToPath.resize(InputSeqCount);

#if DEBUG
	vector<bool> Done(InputSeqCount, false);
#endif
	vector<uint> Order;
	InputSeqs.GetLengthOrder(Order);
	uint LastLength = UINT_MAX;
	const float MinEE = (1 - m_MinEA);
	for (uint k = 0; k < InputSeqCount; ++k)
		{
		uint SeqIndex = Order[k];
		asserta(SeqIndex < InputSeqCount);
#if DEBUG
		asserta(Done[SeqIndex] == false);
		Done[SeqIndex] = true;
#endif
		const uint L = (uint) InputSeqs.GetSeqLength(SeqIndex);
		asserta(L <= LastLength);
		LastLength = L;

		ProgressStep(k, InputSeqCount,
		  "UCLUST %u seqs EE<%.2f, %u centroids, %u members",
		  InputSeqCount, MinEE, CentroidCount, MemberCount);

		string &Path = m_SeqIndexToPath[SeqIndex];
		uint RepSeqIndex = Search(SeqIndex, Path);
		if (RepSeqIndex == UINT_MAX)
			{
			m_CentroidSeqIndexes.push_back(SeqIndex);
			AddSeqToIndex(SeqIndex);
			++CentroidCount;
			RepSeqIndex = SeqIndex;
			Path.clear();
			}
		else
			++MemberCount;
		  
		m_SeqIndexToCentroidSeqIndex[SeqIndex] = RepSeqIndex;
		}
	}

void UClust::GetCentroidSeqs(MultiSequence &CentroidSeqs) const
	{
	const uint CentroidCount = SIZE(m_CentroidSeqIndexes);
	for (uint i = 0; i < CentroidCount; ++i)
		{
		uint SeqIndex = m_CentroidSeqIndexes[i];
		const Sequence *Seq = m_InputSeqs->GetSequence(SeqIndex);
		CentroidSeqs.AddSequence(Seq, false);
		}
	AssertSameLabels(CentroidSeqs);
	}

void UClust::GetGSIs(
  vector<uint> &CentroidGSIs,
  vector<uint> &MemberGSIs,
  vector<uint> &MemberCentroidGSIs,
  vector<string> &GSIToMemberCentroidPath) const
	{
	CentroidGSIs.clear();
	MemberGSIs.clear();
	MemberCentroidGSIs.clear();
	GSIToMemberCentroidPath.clear();

	const uint InputSeqCount = m_InputSeqs->GetSeqCount();
	const uint GSICount = GetGSICount();

	GSIToMemberCentroidPath.resize(GSICount);

	const uint ClusterCount = SIZE(m_CentroidSeqIndexes);
	for (uint ClusterIndex = 0; ClusterIndex < ClusterCount; ++ClusterIndex)
		{
		uint CentroidSeqIndex = m_CentroidSeqIndexes[ClusterIndex];
		const Sequence *Seq = m_InputSeqs->GetSequence(CentroidSeqIndex);
		uint CentroidGSI = GetGSIByLabel(Seq->m_Label);
		CentroidGSIs.push_back(CentroidGSI);
		}

	asserta(SIZE(m_SeqIndexToCentroidSeqIndex) == InputSeqCount);
	for (uint MemberSeqIndex = 0; MemberSeqIndex < InputSeqCount; ++MemberSeqIndex)
		{
		uint CentroidSeqIndex = m_SeqIndexToCentroidSeqIndex[MemberSeqIndex];
		if (CentroidSeqIndex == MemberSeqIndex)
			continue;

		const string &Path = m_SeqIndexToPath[MemberSeqIndex];
		const Sequence *MemberSeq = m_InputSeqs->GetSequence(MemberSeqIndex);
		const Sequence *CentroidSeq = m_InputSeqs->GetSequence(CentroidSeqIndex);

		uint MemberGSI = GetGSIByLabel(MemberSeq->m_Label);
		uint MemberCentroidGSI = GetGSIByLabel(CentroidSeq->m_Label);

		MemberGSIs.push_back(MemberGSI);
		MemberCentroidGSIs.push_back(MemberCentroidGSI);

		asserta(!Path.empty());
		GSIToMemberCentroidPath[MemberGSI] = Path;
		}
	}

void cmd_uclust()
	{
	const string &InputFileName = opt(uclust);
	const string &OutputFileName = opt(output);
	const float MinEA = (float) optd(minea, 0.9);

	//LoadGlobalInputMS(InputFileName);
	//MultiSequence &InputSeqs = GetGlobalInputMS();
	MultiSequence InputSeqs;
	LoadInput(InputSeqs);

	bool IsNucleo = InputSeqs.GuessIsNucleo();
	if (IsNucleo)
		SetAlpha(ALPHA_Nucleo);
	else
		SetAlpha(ALPHA_Amino);

	UClust U;
	U.Run(InputSeqs, MinEA);

	MultiSequence *CentroidSeqs = new MultiSequence;
	U.GetCentroidSeqs(*CentroidSeqs);
	CentroidSeqs->WriteMFA(OutputFileName);
	}
