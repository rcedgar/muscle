#include "muscle.h"
#include "mpcflat.h"

void MPCFlat::RefineIter()
	{
	set<int> SeqIndexes1, SeqIndexes2;

	const uint SeqCount = GetSeqCount();
	asserta(m_MSA != 0);
	asserta(m_MSA->GetSeqCount() == SeqCount);

	// create two separate groups
	for (uint SeqIndex = 0; SeqIndex < SeqCount; SeqIndex++)
		if (rand()%2 == 0)
			SeqIndexes1.insert(SeqIndex);
		else
			SeqIndexes2.insert(SeqIndex);

	if (SeqIndexes1.empty() || SeqIndexes2.empty())
		return;

	const MultiSequence *MSA1 = m_MSA->Project(SeqIndexes1);
	const MultiSequence *MSA2 = m_MSA->Project(SeqIndexes2);
	delete m_MSA;

	MultiSequence *MSA12 = AlignAlns(*MSA1, *MSA2);
	m_MSA = MSA12;

	delete MSA1;
	delete MSA2;
	}
