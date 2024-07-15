#include "muscle.h"
#include "tree.h"
#include "mpcflat.h"

void MPCFlat::FreeProgMSAs()
	{
	const uint n = SIZE(m_ProgMSAs);
	for (uint i = 0; i < n; ++i)
		{
		MultiSequence *MSA = m_ProgMSAs[i];
		if (MSA != 0)
			delete MSA;
		}
	m_ProgMSAs.clear();
	}

void MPCFlat::FreeSparsePosts()
	{
	for (uint i = 0; i < SIZE(m_SparsePosts1); ++i)
		{
		if (m_SparsePosts1[i] != 0)
			{
			delete m_SparsePosts1[i];
			m_SparsePosts1[i] = 0;
			}
		}

	for (uint i = 0; i < SIZE(m_SparsePosts2); ++i)
		{
		if (m_SparsePosts2[i] != 0)
			{
			delete m_SparsePosts2[i];
			m_SparsePosts2[i] = 0;
			}
		}

	m_SparsePosts1.clear();
	m_SparsePosts2.clear();
	}

void MPCFlat::ProgAln(uint JoinIndex)
	{
	uint Index1 = m_JoinIndexes1[JoinIndex];
	uint Index2 = m_JoinIndexes2[JoinIndex];
	assert(Index1 < SIZE(m_ProgMSAs));
	assert(Index2 < SIZE(m_ProgMSAs));

	MultiSequence *MSA1 = m_ProgMSAs[Index1];
	MultiSequence *MSA2 = m_ProgMSAs[Index2];
	assert(MSA1 != 0);
	assert(MSA2 != 0);
	MultiSequence *MSA12 = AlignAlns(*MSA1, *MSA2);

#if 0//TRACE
	uint SeqCount12 = MSA12->GetSeqCount();
	uint Index12 = SIZE(m_ProgMSAs);

	uint SeqCount1 = MSA1->GetSeqCount();
	uint SeqCount2 = MSA2->GetSeqCount();
	Log("Flat Join %u(%u) + %u(%u) = %u(%u)\n",
	  Index1, SeqCount1, Index2, SeqCount2, Index12, SeqCount12);
#endif

	m_ProgMSAs.push_back(MSA12);
	const uint SeqCount = GetSeqCount();
	delete MSA1;
	delete MSA2;

	m_ProgMSAs[Index1] = 0;
	m_ProgMSAs[Index2] = 0;
	}

void MPCFlat::ProgressiveAlign()
	{
	const uint SeqCount = m_MyInputSeqs->GetSeqCount();
	const uint JoinCount = SeqCount - 1;
	const uint NodeCount = SeqCount + JoinCount;

	for (uint i = 0; i < SeqCount; ++i)
		{
		const Sequence *Seq = m_MyInputSeqs->GetSequence(i);
		MultiSequence *MS = new MultiSequence;
		MS->AddSequence(Seq, false);
		m_ProgMSAs.push_back(MS);
		}

	asserta(SIZE(m_JoinIndexes1) == JoinCount);
	asserta(SIZE(m_JoinIndexes2) == JoinCount);

	ValidateJoinOrder(m_JoinIndexes1, m_JoinIndexes2);

	for (uint JoinIndex = 0; JoinIndex < JoinCount; ++JoinIndex)
		ProgAln(JoinIndex);

	asserta(SIZE(m_ProgMSAs) == NodeCount);
	m_MSA = m_ProgMSAs[NodeCount-1];
	m_ProgMSAs[NodeCount-1] = 0;
	FreeProgMSAs();
	asserta(m_MSA != 0);
	}
