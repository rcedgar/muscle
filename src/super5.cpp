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

void Super5::SetOpts()
	{
	m_MinEAPass1 = (float) optd(super5_minea1, DEFAULT_MIN_EA_SUPER5_PASS1);
	}

void Super5::ClearTreesAndMSAs()
	{
	m_GuideTree_None.Clear();
	m_GuideTree_ABC.Clear();
	m_GuideTree_ACB.Clear();
	m_GuideTree_BCA.Clear();

	m_FinalMSA_None.Clear();
	m_FinalMSA_ABC.Clear();
	m_FinalMSA_ACB.Clear();
	m_FinalMSA_BCA.Clear();
	}

void Super5::MakeCentroidSeqs(MultiSequence &InputSeqs)
	{
	m_InputSeqs = &InputSeqs;
	m_UniqueSeqs = new MultiSequence;
	m_CentroidSeqs = new MultiSequence;
	m_CentroidMSA = new MultiSequence;

	m_D.Run(*m_InputSeqs);
	m_D.Validate();
	m_D.GetUniqueSeqs(*m_UniqueSeqs);
	SetDupeVecs();

	m_U.Run(*m_UniqueSeqs, m_MinEAPass1);
	m_U.GetCentroidSeqs(*m_CentroidSeqs);
	SetCentroidVecs();
	SetCentroidSeqsVecs();
	ValidateVecs();
	}

void Super5::Run(MultiSequence &InputSeqs, TREEPERM Perm)
	{
	MakeCentroidSeqs(InputSeqs);
	m_S4.Run(*m_CentroidSeqs, Perm);

	if (Perm != TP_All)
		{
		m_CentroidMSA = &m_S4.m_FinalMSA;
		SetCentroidMSAVecs();
		AlignMembers();
		AlignDupes();
		m_FinalMSA = m_ExtendedMSA;
		SortMSA(*m_FinalMSA);
		return;
		}

	m_CentroidMSA = &m_S4.m_FinalMSA_None;
	SetCentroidMSAVecs();
	AlignMembers();
	AlignDupes();
	m_FinalMSA_None.Copy(*m_ExtendedMSA);
	SortMSA(m_FinalMSA_None);

	m_CentroidMSA = &m_S4.m_FinalMSA_ABC;
	SetCentroidMSAVecs();
	AlignMembers();
	AlignDupes();
	m_FinalMSA_ABC.Copy(*m_ExtendedMSA);
	SortMSA(m_FinalMSA_ABC);

	m_CentroidMSA = &m_S4.m_FinalMSA_ACB;
	SetCentroidMSAVecs();
	AlignMembers();
	AlignDupes();
	m_FinalMSA_ACB.Copy(*m_ExtendedMSA);
	SortMSA(m_FinalMSA_ACB);

	m_CentroidMSA = &m_S4.m_FinalMSA_BCA;
	SetCentroidMSAVecs();
	AlignMembers();
	AlignDupes();
	m_FinalMSA_BCA.Copy(*m_ExtendedMSA);
	SortMSA(m_FinalMSA_BCA);
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
		uint GSI = GetGSIByLabel(Seq->m_Label);
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
		uint GSI = GetGSIByLabel(Seq->m_Label);
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
	m_DupeRepGSIToMemberGSIs.clear();
	m_DupeRepGSIToMemberGSIs.resize(InputSeqCount);

	m_D.GetDupeGSIs(
	  m_DupeGSIs,
	  m_DupeRepGSIs);

	m_IsDupe.resize(InputSeqCount, false);
	const uint DupeCount = SIZE(m_DupeGSIs);
	for (uint i = 0; i < DupeCount; ++i)
		{
		uint GSI = m_DupeGSIs[i];
		uint DupeRepGSI = m_DupeRepGSIs[i];
		asserta(GSI < InputSeqCount);
		asserta(m_IsDupe[GSI] == false);
		m_IsDupe[GSI] = true;
		m_DupeRepGSIToMemberGSIs[DupeRepGSI].push_back(GSI);
		}
	}

void Super5::SetCentroidVecs()
	{
	const uint InputSeqCount = m_InputSeqs->GetSeqCount();

	m_CentroidGSIs.clear();
	m_MemberGSIs.clear();
	m_MemberCentroidGSIs.clear();

	m_GSIToCentroidGSI.clear();
	m_GSIToCentroidGSI.resize(InputSeqCount, UINT_MAX);

	m_CentroidGSIToMemberGSIs.clear();
	m_CentroidGSIToMemberGSIs.resize(InputSeqCount);

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
		m_GSIToCentroidGSI[MemberGSI] = MemberCentroidGSI;
		m_CentroidGSIToMemberGSIs[MemberCentroidGSI].push_back(MemberGSI);
		m_GSIToMemberCount[MemberCentroidGSI] += 1;
		}
	}

void Super5::LogClusters() const
	{
	const uint InputSeqCount = m_InputSeqs->GetSeqCount();
	asserta(SIZE(m_IsDupe) == InputSeqCount);
	asserta(SIZE(m_IsCentroid) == InputSeqCount);
	asserta(SIZE(m_IsMember) == InputSeqCount);
	asserta(SIZE(m_GSIToCentroidSeqsSeqIndex) == InputSeqCount);
	asserta(SIZE(m_GSIToCentroidMSASeqIndex) == InputSeqCount);
	asserta(SIZE(m_GSIToMemberCount) == InputSeqCount);

	Log("%5.5s", "GSI");
	Log("  %3.3s", "Cat");
	Log("  %5.5s", "CSSI");
	Log("  %5.5s", "CMSI");
	Log("  %5.5s", "MbCt");
	Log("  %5.5s", "Cent");
	Log("\n");

	for (uint GSI = 0; GSI < InputSeqCount; ++GSI)
		{
		string Label;
		GetLabelByGSI(GSI, Label);
		bool Dupe = m_IsDupe[GSI];
		bool Centroid = m_IsCentroid[GSI];
		bool Member = m_IsMember[GSI];
		const char *Cat = "Error";
		if (Dupe)
			Cat = "Dup";
		else if (Centroid)
			Cat = "Cnt";
		else if (Member)
			Cat = "Mem";
		else
			asserta(false);

		uint CentroidSeqsSeqIndex = m_GSIToCentroidSeqsSeqIndex[GSI];
		uint CentroidMSASeqIndex = m_GSIToCentroidMSASeqIndex[GSI];
		uint MemberCount = m_GSIToMemberCount[GSI];
		uint CentroidGSI = m_GSIToCentroidGSI[GSI];
		string CentroidLabel = "";
		if (CentroidGSI != UINT_MAX)
			GetLabelByGSI(CentroidGSI, CentroidLabel);

		Log("%5u", GSI);
		Log("  %3.3s", Cat);

		if (CentroidSeqsSeqIndex == UINT_MAX)
			Log("  %5.5s", "*");
		else
			Log("  %5u", CentroidSeqsSeqIndex);

		if (CentroidMSASeqIndex == UINT_MAX)
			Log("  %5.5s", "*");
		else
			Log("  %5u", CentroidMSASeqIndex);

		Log("  %5u", MemberCount);

		if (CentroidGSI == UINT_MAX)
			Log("  %5.5s", "*");
		else
			Log("  %5u", CentroidGSI);

		Log("  >%s", Label.c_str());
		if (CentroidLabel != "")
			Log("  >> %s", CentroidLabel.c_str());

		const vector<uint> &Members = m_CentroidGSIToMemberGSIs[GSI];
		const uint M = SIZE(Members);
		for (uint m = 0; m < min(M, 4u); ++m)
			Log(" <%u>", Members[m]);
		if (M > 4)
			Log(" ...");

		const vector<uint> &DupeMembers = m_DupeRepGSIToMemberGSIs[GSI];
		const uint DM = SIZE(DupeMembers);
		for (uint m = 0; m < min(DM, 4u); ++m)
			Log(" =%u", DupeMembers[m]);
		if (DM > 4)
			Log(" ...");

		Log("\n");
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
		Sequence *MemberSeq = (Sequence *) &GetGlobalInputSeqByIndex(MemberGSI);
		MemberSeqs->AddSequence(MemberSeq, false);

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
		uint GSI = GetGSIByLabel(Seq->m_Label);
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
		string Label;
		GetLabelByGSI(DupeGSI, Label);
		AlignedDupe->OverwriteLabel(Label);
		m_ExtendedMSA->AddSequence(AlignedDupe, true);
		}
	ProgressLog(" done.\n");
	AssertSeqsEqInput(*m_ExtendedMSA);
	}

void Super5::GetLabelToAlnSeqIndex(const MultiSequence &Aln,
  map<string, uint> &LabelToAlnSeqIndex) const
	{
	LabelToAlnSeqIndex.clear();
	const uint SeqCount = Aln.GetSeqCount();
	for (uint i = 0; i < SeqCount; ++i)
		{
		const string &Label = Aln.GetLabelStr(i);
		if (LabelToAlnSeqIndex.find(Label) != LabelToAlnSeqIndex.end())
			Die("Duplicate label >%s", Label.c_str());
		LabelToAlnSeqIndex[Label] = i;
		}
	}

void Super5::SortMSA(MultiSequence &Aln)
	{
	if (opt(input_order))
		SortMSA_ByInputOrder(Aln);
	else
		{
		opt(tree_order);
		SortMSA_ByGuideTree(Aln);
		}
	}

const string &Super5::GetLabel(uint GSI) const
	{
	const string &Label = m_InputSeqs->GetLabelStr(GSI);
	return Label;
	}

void Super5::AppendLabels(uint GSI, vector<string> &Labels) const
	{
	const string &Label = GetLabel(GSI);
	Labels.push_back(Label);
	asserta(GSI < SIZE(m_DupeRepGSIToMemberGSIs));
	const vector<uint> &DupeGSIs = m_DupeRepGSIToMemberGSIs[GSI];
	const uint D = SIZE(DupeGSIs);
	for (uint d = 0; d < D; ++d)
		{
		const string &Label = GetLabel(DupeGSIs[d]);
		Labels.push_back(Label);
		}
	}

void Super5::AppendLabelsFromCentroid(uint CentroidIndex,
  vector<string> &Labels) const
	{
	asserta(CentroidIndex < SIZE(m_CentroidGSIs));
	uint CentroidGSI = m_CentroidGSIs[CentroidIndex];
	AppendLabels(CentroidGSI, Labels);
	asserta(CentroidGSI < SIZE(m_CentroidGSIToMemberGSIs));
	const vector<uint> &MemberGSIs = m_CentroidGSIToMemberGSIs[CentroidGSI];
	const uint M = SIZE(MemberGSIs);
	for (uint m = 0; m < M; ++m)
		AppendLabels(MemberGSIs[m], Labels);
	}

void Super5::GetLabelsInGuideTreeOrder(vector<string> &Labels) const
	{
	Labels.clear();
	const uint CentroidCount = SIZE(m_CentroidGSIs);
	for (uint CentroidIndex = 0; CentroidIndex < CentroidCount; ++CentroidIndex)
		AppendLabelsFromCentroid(CentroidIndex, Labels);

	const uint SeqCount = m_InputSeqs->GetSeqCount();
	asserta(SIZE(Labels) == SeqCount);
	}

void Super5::SortMSA_ByGuideTree(MultiSequence &Aln)
	{
	const uint SeqCount = m_InputSeqs->GetSeqCount();

	map<string, uint> LabelToMSASeqIndex;
	GetLabelToAlnSeqIndex(Aln, LabelToMSASeqIndex);

	vector<string> Labels;
	GetLabelsInGuideTreeOrder(Labels);
	asserta(SIZE(Labels) == SeqCount);

	vector<const Sequence *> SortedSeqs;
	for (uint i = 0; i < SeqCount; ++i)
		{
		const string &Label = Labels[i];
		map<string, uint>::const_iterator p =
		  LabelToMSASeqIndex.find(Label);
		if (p == LabelToMSASeqIndex.end())
			Die("Super5::SortMSA_ByGuideTree(), missing >%s", Label.c_str());
		uint MSASeqIndex = p->second;
		const Sequence *Seq = Aln.GetSequence(MSASeqIndex);
		SortedSeqs.push_back(Seq);
		}

	for (uint i = 0; i < SeqCount; ++i)
		Aln.m_Seqs[i] = SortedSeqs[i];
	}

void Super5::SortMSA_ByInputOrder(MultiSequence &Aln)
	{
	const uint SeqCount = m_InputSeqs->GetSeqCount();
	map<string, uint> LabelToMSASeqIndex;
	GetLabelToAlnSeqIndex(Aln, LabelToMSASeqIndex);

	vector<const Sequence *> SortedSeqs;
	for (uint i = 0; i < SeqCount; ++i)
		{
		const string &Label = m_InputSeqs->GetLabelStr(i);
		map<string, uint>::const_iterator p =
		  LabelToMSASeqIndex.find(Label);
		if (p == LabelToMSASeqIndex.end())
			Die("Super5::SortMSA_ByInputOrder(), missing >%s", Label.c_str());
		uint MSASeqIndex = p->second;
		const Sequence *Seq = Aln.GetSequence(MSASeqIndex);
		SortedSeqs.push_back(Seq);
		}

	for (uint i = 0; i < SeqCount; ++i)
		Aln.m_Seqs[i] = SortedSeqs[i];
	}

void cmd_super5()
	{
	if (optset_mega)
		Die("-super5 does not support -mega, use -super7");

	//LoadGlobalInputMS(opt(super5));
	MultiSequence InputSeqs;
	LoadInput(InputSeqs);

	string &OutputPattern = opt(output);
	if (OutputPattern.empty())
		Die("Must set -output");

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

	if (optset_diversified)
		Die("-diversified not supported");
	if (optset_replicates)
		Die("-replicates not supported");
	if (optset_stratified)
		Die("-stratified not supported");

	TREEPERM Perm = TP_None;
	if (optset_perm)
		Perm = StrToTREEPERM(opt(perm));

	if (Perm == TP_All && OutputPattern.find('@') == string::npos)
		Die("Must be '@' in output filename with -perm all");

	Super5 S5;
	S5.SetOpts();
	S5.Run(InputSeqs, Perm);
	if (Perm == TP_All)
		{
		string FileName_None;
		string FileName_ABC;
		string FileName_ACB;
		string FileName_BCA;

		uint PerturbSeed = 0;
		if (optset_perturb)
			PerturbSeed = opt(perturb);
		MakeReplicateFileName(OutputPattern, TP_None, PerturbSeed, FileName_None);
		MakeReplicateFileName(OutputPattern, TP_ABC, PerturbSeed, FileName_ABC);
		MakeReplicateFileName(OutputPattern, TP_ACB, PerturbSeed, FileName_ACB);
		MakeReplicateFileName(OutputPattern, TP_BCA, PerturbSeed, FileName_BCA);

		S5.m_FinalMSA_None.WriteMFA(FileName_None);
		S5.m_FinalMSA_ABC.WriteMFA(FileName_ABC);
		S5.m_FinalMSA_ACB.WriteMFA(FileName_ACB);
		S5.m_FinalMSA_BCA.WriteMFA(FileName_BCA);
		}
	else
		{
		uint PerturbSeed = 0;
		if (optset_perturb)
			PerturbSeed = opt(perturb);
		
		string OutputFileName;
		if (OutputPattern.find('@') == string::npos)
			OutputFileName = OutputPattern;
		else
			MakeReplicateFileName(OutputPattern, Perm, PerturbSeed, OutputFileName);

		S5.m_FinalMSA->WriteMFA(OutputFileName);
		}

	if (S5.m_FinalMSA != 0)
		{
		S5.m_FinalMSA->Clear();
		S5.m_FinalMSA = 0;
		}
	S5.m_FinalMSA_None.Clear();
	S5.m_FinalMSA_ABC.Clear();
	S5.m_FinalMSA_ACB.Clear();
	S5.m_FinalMSA_BCA.Clear();
#if SEQ_TRACE
	Sequence::AllocReport("final");
#endif
	}
