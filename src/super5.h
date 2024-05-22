#pragma once

#include "derep.h"
#include "uclust.h"
#include "transaln.h"
#include "super4.h"

static const float DEFAULT_MIN_EA_SUPER5_PASS1 = 0.99f;

class Super5
	{
public:
	float m_MinEAPass1 = DEFAULT_MIN_EA_SUPER5_PASS1;
	MultiSequence *m_InputSeqs = 0;
	MultiSequence *m_UniqueSeqs = 0;
	MultiSequence *m_CentroidSeqs = 0;
	MultiSequence *m_CentroidMSA = 0;
	MultiSequence *m_ExtendedMSA = 0;
	MultiSequence *m_FinalMSA = 0;

	Tree m_GuideTree_None;
	Tree m_GuideTree_ABC;
	Tree m_GuideTree_ACB;
	Tree m_GuideTree_BCA;

	MultiSequence m_FinalMSA_None;
	MultiSequence m_FinalMSA_ABC;
	MultiSequence m_FinalMSA_ACB;
	MultiSequence m_FinalMSA_BCA;

	Derep m_D;
	UClust m_U;
	TransAln m_TA;
	Super4 m_S4;

	vector<bool> m_IsDupe;
	vector<bool> m_IsCentroid;
	vector<bool> m_IsMember;

	vector<uint> m_DupeGSIs;
	vector<uint> m_DupeRepGSIs;

	vector<uint> m_CentroidGSIs;
	vector<uint> m_MemberGSIs;
	vector<uint> m_MemberCentroidGSIs;

	vector<uint> m_CentroidSeqsSeqIndexToGSI;
	vector<uint> m_CentroidMSASeqIndexToGSI;

	vector<uint> m_GSIToCentroidSeqsSeqIndex;
	vector<uint> m_GSIToCentroidMSASeqIndex;
	vector<uint> m_GSIToMemberCount;
	vector<uint> m_GSIToCentroidGSI;
	vector<vector<uint> > m_CentroidGSIToMemberGSIs;
	vector<vector<uint> > m_DupeRepGSIToMemberGSIs;

	vector<string> m_GSIToMemberCentroidPath;

public:
	void SetOpts();
	void Run(MultiSequence &InputSeqs, TREEPERM Perm);
	void MakeCentroidSeqs(MultiSequence &InputSeqs);
	void SetDupeVecs();
	void SetCentroidVecs();
	void SetCentroidSeqsVecs();
	void SetCentroidMSAVecs();
	void AlignMembers();
	void AlignDupes();
	void ValidateVecs() const;
	void ClearTreesAndMSAs();
	void LogClusters() const;

	void GetLabelsInGuideTreeOrder(vector<string> &Labels) const;
	void AppendLabelsFromCentroid(uint CentroidIndex,
	  vector<string> &Labels) const;
	const string &GetLabel(uint GSI) const;
	void AppendLabels(uint GSI, vector<string> &Labels) const;
	void GetLabelToAlnSeqIndex(const MultiSequence &Aln,
	  map<string, uint> &LabelToAlnSeqIndex) const;
	void SortMSA_ByInputOrder(MultiSequence &Aln);
	void SortMSA_ByGuideTree(MultiSequence &Aln);
	void SortMSA(MultiSequence &Aln);
	};
