#pragma once

#include "derep.h"
#include "uclust.h"
#include "transaln.h"

class Super5
	{
public:
	float m_UClustEA = 0.99f;
	MultiSequence *m_InputSeqs = 0;
	MultiSequence *m_UniqueSeqs = 0;
	MultiSequence *m_CentroidSeqs = 0;
	MultiSequence *m_CentroidMSA = 0;
	MultiSequence *m_ExtendedMSA = 0;
	MultiSequence *m_FinalMSA = 0;
	TREEPERM m_TreePerm = TP_None;

	Derep m_D;
	UClust m_U;
	TransAln m_TA;

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

	vector<string> m_GSIToMemberCentroidPath;

public:
	void Run(MultiSequence &InputSeqs);
	void SetDupeVecs();
	void SetCentroidVecs();
	void SetCentroidSeqsVecs();
	void SetCentroidMSAVecs();
	void AlignMembers();
	void AlignDupes();
	void ValidateVecs() const;
	};
