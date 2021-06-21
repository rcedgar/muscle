#include "muscle.h"
#include "tree.h"
#include "mpc.h"

void MPC::ProgressiveAlign()
	{
	vector<MultiSequence *> MSAs;
	const uint SeqCount = m_InputSeqs->GetSeqCount();
	const uint JoinCount = SeqCount - 1;
	const uint NodeCount = SeqCount + JoinCount;

	for (uint i = 0; i < SeqCount; ++i)
		{
		Sequence *Seq = m_InputSeqs->GetSequence(i);
		MultiSequence *MS = new MultiSequence;
		MS->AddSequence(Seq);
		MSAs.push_back(MS);
		}

	asserta(SIZE(m_JoinIndexes1) == JoinCount);
	asserta(SIZE(m_JoinIndexes2) == JoinCount);

	ValidateJoinOrder(m_JoinIndexes1, m_JoinIndexes2);

	for (uint JoinIndex = 0; JoinIndex < JoinCount; ++JoinIndex)
		{
		uint Index1 = m_JoinIndexes1[JoinIndex];
		uint Index2 = m_JoinIndexes2[JoinIndex];
		asserta(Index1 < SIZE(MSAs));
		asserta(Index2 < SIZE(MSAs));

		MultiSequence *MSA1 = MSAs[Index1];
		MultiSequence *MSA2 = MSAs[Index2];
		asserta(MSA1 != 0);
		asserta(MSA2 != 0);

		MultiSequence *MSA12 = AlignAlignments(MSA1, MSA2, m_SparseMatrices);
		MSAs.push_back(MSA12);

		if (Index1 >= SeqCount)
			{
			delete MSA1;
			MSAs[Index1] = 0;
			}
		if (Index2 >= SeqCount)
			{
			delete MSA2;
			MSAs[Index2] = 0;
			}
		}
	asserta(SIZE(MSAs) == NodeCount);
	m_MSA = MSAs[NodeCount-1];
	asserta(m_MSA != 0);
	}
