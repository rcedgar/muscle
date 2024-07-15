#include "muscle.h"
#include "derep.h"

void Derep::Clear()
	{
	m_SeqIndexToRepSeqIndex.clear();
	m_RepSeqIndexes.clear();
	m_RepSeqIndexToSeqIndexes.clear();
	m_HashToSeqIndexes.clear();
	}

// FNV64 hash
uint Derep::CalcHash(const Sequence *Seq) const
	{
    uint64 hash = 0xcbf29ce484222325uL;
	const uint L = Seq->GetLength();
    for (uint i = 0; i < L; ++i)
	    {
		char c = Seq->GetChar(i);
        byte b = (byte) tolower(c);
        hash *= 1099511628211uL;
        hash ^= (uint64) b;
		}
	uint h = uint(hash%m_SlotCount);
	return h;
	}

void Derep::Run(MultiSequence &InputSeqs, bool ShowProgress)
	{
	Clear();
	m_InputSeqs = &InputSeqs;
	const uint InputSeqCount = InputSeqs.GetSeqCount();
	m_SlotCount = 3*InputSeqCount + 7;

	m_HashToSeqIndexes.resize(m_SlotCount);
	m_SeqIndexToRepSeqIndex.resize(InputSeqCount, UINT_MAX);
	m_RepSeqIndexToSeqIndexes.resize(InputSeqCount);

	uint UniqueCount = 0;
	uint DupeCount = 0;
	for (uint SeqIndex = 0; SeqIndex < InputSeqCount; ++SeqIndex)
		{
		uint RepSeqIndex = Search(SeqIndex);
		if (RepSeqIndex == UINT_MAX)
			{
			AddToHash(SeqIndex);
			asserta(SIZE(m_RepSeqIndexes) == UniqueCount);
			m_RepSeqIndexes.push_back(SeqIndex);
			m_RepSeqIndexToSeqIndexes[SeqIndex].push_back(SeqIndex);
			m_SeqIndexToRepSeqIndex[SeqIndex] = SeqIndex;
			++UniqueCount;
			}
		else
			{
			m_RepSeqIndexToSeqIndexes[RepSeqIndex].push_back(SeqIndex);
			m_SeqIndexToRepSeqIndex[SeqIndex] = RepSeqIndex;
			++DupeCount;
			}
		if (ShowProgress)
			ProgressStep(SeqIndex, InputSeqCount, "Derep %u uniques, %u dupes",
				UniqueCount, DupeCount);
		}
	}

bool Derep::SeqsEq(uint SeqIndex1, uint SeqIndex2) const
	{
	if (m_Disable)
		return false;
	const Sequence *Seq1 = m_InputSeqs->GetSequence(SeqIndex1);
	const Sequence *Seq2 = m_InputSeqs->GetSequence(SeqIndex2);
	const uint L = Seq1->GetLength();
	const uint L2 = Seq2->GetLength();
	if (L2 != L)
		return false;
	for (uint i = 0; i < L; ++i)
		{
		char c1 = Seq1->GetChar(i);
		char c2 = Seq2->GetChar(i);
		if (toupper(c1) != toupper(c2))
			return false;
		}
	return true;
	}

uint Derep::Search(uint SeqIndex) const
	{
	if (m_Disable)
		return UINT_MAX;
	const Sequence *Seq = m_InputSeqs->GetSequence(SeqIndex);
	asserta(Seq != 0);
	uint h = CalcHash(Seq);
	asserta(h < SIZE(m_HashToSeqIndexes));
	const vector<uint> &Row = m_HashToSeqIndexes[h];
	const uint n = SIZE(Row);
	for (uint i = 0; i < n; ++i)
		{
		uint SeqIndex2 = Row[i];
		if (SeqsEq(SeqIndex, SeqIndex2))
			return SeqIndex2;
		}
	return UINT_MAX;
	}

void Derep::AddToHash(uint SeqIndex)
	{
	const Sequence *Seq = m_InputSeqs->GetSequence(SeqIndex);
	asserta(Seq != 0);
	uint h = CalcHash(Seq);
	asserta(h < SIZE(m_HashToSeqIndexes));
	vector<uint> &Row = m_HashToSeqIndexes[h];
	Row.push_back(SeqIndex);
	}

void Derep::GetUniqueSeqs(MultiSequence &UniqueSeqs)
	{
	asserta(UniqueSeqs.GetSeqCount() == 0);
	const uint UniqueCount = SIZE(m_RepSeqIndexes);
	for (uint i = 0; i < UniqueCount; ++i)
		{
		uint SeqIndex = m_RepSeqIndexes[i];
		const Sequence *Seq = m_InputSeqs->GetSequence(SeqIndex);
		UniqueSeqs.AddSequence(Seq, false);
		}
	//AssertSameLabels(UniqueSeqs);
	}

void Derep::GetRepLabelToDupeLabels(
  unordered_map<string, vector<string> > &RepLabelToMemberLabels) const
	{
	RepLabelToMemberLabels.clear();
	const uint SeqCount = m_InputSeqs->GetSeqCount();
	const uint ClusterCount = SIZE(m_RepSeqIndexes);

	vector<string> Empty;
	for (uint ClusterIndex = 0; ClusterIndex < ClusterCount; ++ClusterIndex)
		{
		uint RepSeqIndex = m_RepSeqIndexes[ClusterIndex];
		const string &RepLabel = m_InputSeqs->GetLabelStr(RepSeqIndex);
		const vector<uint> &MemberSeqIndexes =
		  m_RepSeqIndexToSeqIndexes[RepSeqIndex];
		const uint MemberCount = SIZE(MemberSeqIndexes);
		if (MemberCount == 1)
			continue;
		RepLabelToMemberLabels[RepLabel] = Empty;
		for (uint MemberIndex = 0; MemberIndex < MemberCount; ++MemberIndex)
			{
			uint MemberSeqIndex = MemberSeqIndexes[MemberIndex];
			const string &MemberLabel =
			  m_InputSeqs->GetLabelStr(MemberSeqIndex);
			if (MemberLabel != RepLabel)
				RepLabelToMemberLabels[RepLabel].push_back(MemberLabel);
			}
		}
	}

void Derep::Validate() const
	{
	asserta(m_InputSeqs != 0);
	const uint InputSeqCount = m_InputSeqs->GetSeqCount();
	asserta(SIZE(m_SeqIndexToRepSeqIndex) == InputSeqCount);
	asserta(SIZE(m_RepSeqIndexToSeqIndexes) == InputSeqCount);

	const uint ClusterCount = SIZE(m_RepSeqIndexes);

	set<uint> RepSeqIndexSet;
	for (uint SeqIndex = 0; SeqIndex < InputSeqCount; ++SeqIndex)
		{
		uint RepSeqIndex = m_SeqIndexToRepSeqIndex[SeqIndex];
		RepSeqIndexSet.insert(RepSeqIndex);
		}

	const uint RepSeqIndexCount = SIZE(m_RepSeqIndexes);
	asserta(SIZE(RepSeqIndexSet) == RepSeqIndexCount);
	for (uint i = 0; i < RepSeqIndexCount; ++i)
		{
		uint RepSeqIndex = m_RepSeqIndexes[i];
		asserta(RepSeqIndexSet.find(RepSeqIndex) != RepSeqIndexSet.end());

		const vector<uint> &MemberSeqIndexes =
		  m_RepSeqIndexToSeqIndexes[RepSeqIndex];

		const uint MemberCount = SIZE(MemberSeqIndexes);
		asserta(MemberCount > 0);
		for (uint j = 0; j < MemberCount; ++j)
			{
			uint MemberSeqIndex = MemberSeqIndexes[j];
			uint MemberRepSeqIndex = m_SeqIndexToRepSeqIndex[MemberSeqIndex];
			asserta(MemberRepSeqIndex == RepSeqIndex);
			}
		}
	}

void Derep::GetDupeGSIs(vector<uint> &GSIs,
  vector<uint> &GlobalRepSeqIndexes) const
	{
	GSIs.clear();
	GlobalRepSeqIndexes.clear();

	const uint InputSeqCount = m_InputSeqs->GetSeqCount();
	const uint GlobalMSSeqCount = GetGlobalMSSeqCount();
	const uint ClusterCount = SIZE(m_RepSeqIndexes);

	for (uint ClusterIndex = 0; ClusterIndex < ClusterCount; ++ClusterIndex)
		{
		uint RepSeqIndex = m_RepSeqIndexes[ClusterIndex];
		asserta(RepSeqIndex < InputSeqCount);
		const vector<uint> &MemberSeqIndexes =
		  m_RepSeqIndexToSeqIndexes[RepSeqIndex];
		const uint MemberCount = SIZE(MemberSeqIndexes);
		const Sequence *Seq = m_InputSeqs->GetSequence(RepSeqIndex);
		uint GlobalRepSeqIndex = GetGSIByLabel(Seq->m_Label);
		asserta(GlobalRepSeqIndex < GlobalMSSeqCount);
		asserta(MemberSeqIndexes[0] == RepSeqIndex);
		for (uint i = 1; i < MemberCount; ++i)
			{
			uint MemberSeqIndex = MemberSeqIndexes[i];
			const Sequence *Seq = m_InputSeqs->GetSequence(MemberSeqIndex);
			uint GlobalMemberSeqIndex = GetGSIByLabel(Seq->m_Label);
			asserta(GlobalMemberSeqIndex < GlobalMSSeqCount);
			GSIs.push_back(GlobalMemberSeqIndex);
			GlobalRepSeqIndexes.push_back(GlobalRepSeqIndex);
			}
		}
	}

void cmd_derep()
	{
	const string &InputFileName = opt(derep);
	const string &OutputFileName = opt(output);

	MultiSequence InputSeqs;
	InputSeqs.FromFASTA(InputFileName);

	Derep D;
	D.Run(InputSeqs);
	D.Validate();

	MultiSequence *UniqueSeqs = new MultiSequence;
	D.GetUniqueSeqs(*UniqueSeqs);
	UniqueSeqs->WriteMFA(OutputFileName);
	}
