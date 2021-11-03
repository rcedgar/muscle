#include "muscle.h"
#include "mpcflat.h"
#include "locallock.h"

void MPCFlat::ConsIter(uint Iter)
	{
	uint PairCount = SIZE(m_Pairs);
	asserta(PairCount > 0);
	unsigned ThreadCount = GetRequestedThreadCount();
	uint PairCounter = 0;
#pragma omp parallel for num_threads(ThreadCount)
	for (int PairIndex = 0; PairIndex < (int) PairCount; ++PairIndex)
		{
		Lock();
		ProgressStep(PairCounter++, PairCount, "Consistency (%u/%u)",
		  Iter+1, m_ConsistencyIterCount);
		Unlock();

		ConsPair(PairIndex);
		}

	swap(m_ptrSparsePosts, m_ptrUpdatedSparsePosts);
	}
