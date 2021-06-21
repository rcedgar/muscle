#include "muscle.h"
#include "super4.h"
#include "sequence.h"
#include "multisequence.h"
#include "usorter.h"
#include "upgma5.h"
#include "pprog.h"
#include "derep.h"
#include "uclust.h"
#include "super5.h"

void CharVecToStr(const vector<char> &Vec, string &Str)
	{
	Str.clear();
	for (uint i = 0; i < SIZE(Vec); ++i)
		Str += Vec[i];
	}

void Super5::Run(MultiSequence &InputSeqs)
	{
	m_InputSeqs = &InputSeqs;
	m_UniqueSeqs = new MultiSequence;
	m_CentroidSeqs = new MultiSequence;
	m_CentroidMSA = new MultiSequence;

	m_D.Run(*m_InputSeqs);
	m_D.Validate();
	m_D.GetUniqueSeqs(*m_UniqueSeqs);
	SetDupeVecs();

	m_U.Run(*m_UniqueSeqs, m_UClustEA);
	m_U.GetCentroidSeqs(*m_CentroidSeqs);
	SetCentroidVecs();
	SetCentroidSeqsVecs();
	ValidateVecs();

	RunSuper4(*m_CentroidSeqs, *m_CentroidMSA, m_TreePerm);
	SetCentroidMSAVecs();
	AlignMembers();
	AlignDupes();
	m_FinalMSA = m_ExtendedMSA;
	}

void Super5::SetCentroidMSAVecs()
	{
	m_CentroidMSASeqIndexToGSI.clear();
	m_GSIToCentroidMSASeqIndex.clear();

	const uint CentroidSeqCount = SIZE(m_CentroidGSIs);
	const uint CentroidMSASeqCount = m_CentroidMSA->GetSeqCount();
	asserta(CentroidSeqCount == CentroidMSASeqCount);

	const uint GlobalSeqCount = GetGlobalMSSeqCount();
	m_GSIToCentroidMSASeqIndex.resize(GlobalSeqCount, UINT_MAX);
	m_CentroidMSASeqIndexToGSI.clear();

	for (uint CentroidMSASeqIndex = 0; CentroidMSASeqIndex < CentroidSeqCount;
	  ++CentroidMSASeqIndex)
		{
		const Sequence *Seq = m_CentroidMSA->GetSequence(CentroidMSASeqIndex);
		uint GSI = Seq->GetGSI();
		asserta(GSI < GlobalSeqCount);
		if (m_GSIToCentroidMSASeqIndex[GSI] != UINT_MAX)
			Die("Super5::SetCentroidMSAVecs() GSI=%u found twice (%u,%u)",
			  GSI, m_GSIToCentroidMSASeqIndex[GSI],
			  CentroidMSASeqIndex);
		m_GSIToCentroidMSASeqIndex[GSI] = CentroidMSASeqIndex;
		m_CentroidMSASeqIndexToGSI.push_back(GSI);
		}
	}

void Super5::SetCentroidSeqsVecs()
	{
	m_CentroidSeqsSeqIndexToGSI.clear();
	m_GSIToCentroidSeqsSeqIndex.clear();

	const uint CentroidSeqCount = SIZE(m_CentroidGSIs);
	const uint CentroidSeqSeqCount = m_CentroidSeqs->GetSeqCount();
	asserta(CentroidSeqCount == CentroidSeqSeqCount);

	const uint GlobalSeqCount = GetGlobalMSSeqCount();
	m_GSIToCentroidSeqsSeqIndex.resize(GlobalSeqCount, UINT_MAX);
	m_CentroidSeqsSeqIndexToGSI.clear();

	for (uint CentroidSeqSeqIndex = 0; CentroidSeqSeqIndex < CentroidSeqCount;
	  ++CentroidSeqSeqIndex)
		{
		const Sequence *Seq = m_CentroidSeqs->GetSequence(CentroidSeqSeqIndex);
		uint GSI = Seq->GetGSI();
		asserta(GSI < GlobalSeqCount);
		asserta(m_GSIToCentroidSeqsSeqIndex[GSI] == UINT_MAX);
		m_GSIToCentroidSeqsSeqIndex[GSI] = CentroidSeqSeqIndex;
		m_CentroidSeqsSeqIndexToGSI.push_back(GSI);
		}
	}

void Super5::SetDupeVecs()
	{
	const uint InputSeqCount = m_InputSeqs->GetSeqCount();
	m_DupeGSIs.clear();
	m_DupeRepGSIs.clear();
	m_IsDupe.clear();

	m_D.GetDupeGSIs(
	  m_DupeGSIs,
	  m_DupeRepGSIs);

	m_IsDupe.resize(InputSeqCount, false);
	const uint DupeCount = SIZE(m_DupeGSIs);
	for (uint i = 0; i < DupeCount; ++i)
		{
		uint GSI = m_DupeGSIs[i];
		asserta(GSI < InputSeqCount);
		asserta(m_IsDupe[GSI] == false);
		m_IsDupe[GSI] = true;
		}
	}

void Super5::SetCentroidVecs()
	{
	const uint InputSeqCount = m_InputSeqs->GetSeqCount();

	m_CentroidGSIs.clear();
	m_MemberGSIs.clear();
	m_MemberCentroidGSIs.clear();
	m_U.GetGSIs(m_CentroidGSIs,
	  m_MemberGSIs, m_MemberCentroidGSIs,
	  m_GSIToMemberCentroidPath);

	m_IsCentroid.clear();
	m_IsMember.clear();

	m_IsCentroid.resize(InputSeqCount, false);
	m_IsMember.resize(InputSeqCount, false);

	const uint GSICount = GetGlobalMSSeqCount();
	m_GSIToMemberCount.resize(GSICount, 0);

	const uint CentroidCount = SIZE(m_CentroidGSIs);
	for (uint i = 0; i < CentroidCount; ++i)
		{
		uint CentroidGSI = m_CentroidGSIs[i];
		asserta(CentroidGSI < InputSeqCount);
		asserta(!m_IsDupe[CentroidGSI]);
		asserta(!m_IsCentroid[CentroidGSI]);
		m_IsCentroid[CentroidGSI] = true;
		}

	const uint MemberCount = SIZE(m_MemberGSIs);
	asserta(SIZE(m_MemberCentroidGSIs) == MemberCount);
	for (uint i = 0; i < MemberCount; ++i)
		{
		uint MemberGSI = m_MemberGSIs[i];
		uint MemberCentroidGSI = m_MemberCentroidGSIs[i];

		asserta(MemberGSI < GSICount);
		asserta(MemberCentroidGSI < GSICount);
		asserta(m_IsCentroid[MemberCentroidGSI]);

		bool IsDupe = m_IsDupe[MemberGSI];
		bool IsMember = m_IsMember[MemberGSI];
		bool IsCentroid = m_IsCentroid[MemberGSI];
		
		if (IsDupe || IsMember || IsCentroid)
			Die("Super5::SetCentroidVecs(), MemberGSI=%u dupe=%c mem=%c cent=%c",
			   MemberGSI, tof(IsDupe), tof(IsMember), tof(IsCentroid));

		asserta(!IsDupe);
		asserta(!IsMember);
		asserta(!IsCentroid);

		m_IsMember[MemberGSI] = true;

		m_GSIToMemberCount[MemberCentroidGSI] += 1;
		}
	}

void Super5::ValidateVecs() const
	{
	const uint InputSeqCount = m_InputSeqs->GetSeqCount();
	asserta(SIZE(m_IsDupe) == InputSeqCount);
	asserta(SIZE(m_IsCentroid) == InputSeqCount);
	asserta(SIZE(m_IsMember) == InputSeqCount);

	for (uint i = 0; i < InputSeqCount; ++i)
		{
		bool Dupe = m_IsDupe[i];
		bool Centroid = m_IsCentroid[i];
		bool Member = m_IsMember[i];
		int Sum = int(Dupe) + int(Centroid) + int(Member);
		if (Sum != 1)
			Die("Input seq %u dupe %c, centroid %c member %c",
			  i, tof(Dupe), tof(Centroid), tof(Member));
		}
	}

void Super5::AlignMembers()
	{
	const uint MemberCount = SIZE(m_MemberGSIs);
	const uint GSICount = GetGSICount();
	asserta(SIZE(m_MemberCentroidGSIs) == MemberCount);
	asserta(m_CentroidMSA != 0);
	const uint CentroidCount = m_CentroidMSA->GetSeqCount();
	const uint CentroidMSAColCount = m_CentroidMSA->GetColCount();

	vector<uint> MemberIndexToCentroidIndex;
	asserta(SIZE(m_MemberCentroidGSIs) == MemberCount);
	asserta(SIZE(m_GSIToCentroidMSASeqIndex) == GSICount);
	asserta(SIZE(m_GSIToMemberCentroidPath) == GSICount);

	MultiSequence *MemberSeqs = new MultiSequence;
	asserta(MemberSeqs != 0);
	vector<string> MemberPaths;
	for (uint MemberIndex = 0; MemberIndex < MemberCount; ++MemberIndex)
		{
		uint MemberGSI = m_MemberGSIs[MemberIndex];
		Sequence *MemberSeq = (Sequence *) &GetGlobalInputSeq(MemberGSI);
		MemberSeqs->AddSequence(MemberSeq);

		uint CentroidGSI = m_MemberCentroidGSIs[MemberIndex];
		asserta(CentroidGSI < GSICount);
		uint CentroidMSASeqIndex = m_GSIToCentroidMSASeqIndex[CentroidGSI];
		asserta(CentroidMSASeqIndex < CentroidCount);
		MemberIndexToCentroidIndex.push_back(CentroidMSASeqIndex);

		const string &Path = m_GSIToMemberCentroidPath[MemberGSI];
		asserta(!Path.empty());
		MemberPaths.push_back(Path);
		}

	m_TA.Init(*m_CentroidMSA, *MemberSeqs,
	  MemberIndexToCentroidIndex, MemberPaths);
	m_TA.MakeExtendedMSA();
	asserta(m_TA.m_ExtendedMSA != 0);
	m_ExtendedMSA = m_TA.m_ExtendedMSA;
	AssertSeqsEqInput(*m_ExtendedMSA);
	}

void Super5::AlignDupes()
	{
	const uint DupeCount = SIZE(m_DupeGSIs);
	asserta(SIZE(m_DupeRepGSIs) == DupeCount);
	if (DupeCount == 0)
		return;

	ProgressLog("Inserting %u dupes...", DupeCount);
	const uint GSICount = GetGSICount();
	vector<uint> GSIToExtendedSeqIndex(GSICount, UINT_MAX);
	asserta(m_ExtendedMSA != 0);
	const uint ExtendedSeqCount = m_ExtendedMSA->GetSeqCount();
	for (uint ExtendedSeqIndex = 0; ExtendedSeqIndex < ExtendedSeqCount;
	  ++ExtendedSeqIndex)
		{
		const Sequence *Seq = m_ExtendedMSA->GetSequence(ExtendedSeqIndex);
		uint GSI = Seq->GetGSI();
		asserta(GSI < GSICount);
		asserta(GSIToExtendedSeqIndex[GSI] == UINT_MAX);
		GSIToExtendedSeqIndex[GSI] = ExtendedSeqIndex;
		}

	for (uint i = 0; i < DupeCount; ++i)
		{
		uint DupeGSI = m_DupeGSIs[i];
		uint RepGSI = m_DupeRepGSIs[i];
		asserta(RepGSI < GSICount);
		uint RepExtendedSeqIndex = GSIToExtendedSeqIndex[RepGSI];
		asserta(RepExtendedSeqIndex < ExtendedSeqCount);
		const Sequence *Rep = m_ExtendedMSA->GetSequence(RepExtendedSeqIndex);
		Sequence *AlignedDupe = Rep->Clone();
		AlignedDupe->OverwriteGSI(DupeGSI);
		const string &Label = GetGlobalInputSeqLabel(DupeGSI);
		AlignedDupe->OverwriteLabel(Label);
		m_ExtendedMSA->AddSequence(AlignedDupe);
		}
	ProgressLog(" done.\n");
	AssertSeqsEqInput(*m_ExtendedMSA);
	}

void cmd_super5()
	{
	const string &InputFileName = opt(super5);
	const string &OutputFileName = opt(output);
	const float UClustEA = (float) optd(minea, 0.99f);

	MultiSequence &InputSeqs = LoadGlobalInputMS(InputFileName);
	const uint InputSeqCount = GetGlobalMSSeqCount();

	bool IsNucleo = InputSeqs.GuessIsNucleo();
	if (IsNucleo)
		Warning("Input may be nucleotide, a.a. required for -super5");

	SetAlpha(ALPHA_Amino);
	InitProbcons();

	TREEPERM TP = TP_None;
	if (optset_perm)
		TP = StrToTREEPERM(opt(perm));
	if (TP == TP_All)
		Die("-perm all not supported, please specify none, abc, acb or bca");

	Super5 S5;
	S5.m_TreePerm = TP;
	S5.m_UClustEA = UClustEA;
	S5.Run(InputSeqs);

	S5.m_FinalMSA->WriteMFA(OutputFileName);
	}
