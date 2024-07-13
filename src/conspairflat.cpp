#include "muscle.h"
#include "mpcflat.h"

#if 0//TRACE
const byte *g_X;
const byte *g_Y;
const byte *g_Z;
#endif

void MPCFlat::ConsPair(uint PairIndex)
	{
	const pair<uint, uint> &Pair = GetPair(PairIndex);

	uint SeqIndexX = Pair.first;
	uint SeqIndexY = Pair.second;

	const MySparseMx &SparsePostXY = GetSparsePost(PairIndex);

	uint LX = GetSeqLength(SeqIndexX);
	uint LY = GetSeqLength(SeqIndexY);

	asserta(SparsePostXY.GetLX() == LX);
	asserta(SparsePostXY.GetLY() == LY);

	float *Post = AllocPost(LX, LY);
	SparsePostXY.ToPost(Post);

// Account for Z=X and Z=Y (hence the factor 2)
	for (uint k = 0; k < LX*LY; ++k)
		Post[k] *= 2;
#if 0//TRACE
	LogFlatMx("ConsPair Z=X Z=Y", Post, LX, LY);
#endif

	const uint SeqCount = GetSeqCount();
	asserta(SeqIndexX < SeqIndexY); // because convention for pairs
	for (uint SeqIndexZ = 0; SeqIndexZ < SeqCount; ++SeqIndexZ)
		{
		if (SeqIndexZ == SeqIndexX || SeqIndexZ == SeqIndexY)
			continue;
		float wZ = m_Weights[SeqIndexZ];
		wZ = 1.0f;

#if 0//TRACE
		g_X = GetSequence(SeqIndexX)->GetBytePtr();
		g_Y = GetSequence(SeqIndexY)->GetBytePtr();
		g_Z = GetSequence(SeqIndexZ)->GetBytePtr();
#endif
		if (SeqIndexZ < SeqIndexX)
			{
			asserta(SeqIndexZ < SeqIndexY); // because SeqIndexX < SeqIndexY

			uint PairIndexZX = GetPairIndex(SeqIndexZ, SeqIndexX);
			uint PairIndexZY = GetPairIndex(SeqIndexZ, SeqIndexY);

			const MySparseMx &ZX = GetSparsePost(PairIndexZX);
			const MySparseMx &ZY = GetSparsePost(PairIndexZY);

			RelaxFlat_ZX_ZY(ZX, ZY, wZ, Post);
#if 0//TRACE
		LogFlatMx("ConsPair after RelaxFlat_ZX_ZY", Post, LX, LY);
#endif
			}
		else if (SeqIndexZ > SeqIndexX && SeqIndexZ < SeqIndexY)
			{
			uint PairIndexXZ = GetPairIndex(SeqIndexX, SeqIndexZ);
			uint PairIndexZY = GetPairIndex(SeqIndexZ, SeqIndexY);

			const MySparseMx &XZ = GetSparsePost(PairIndexXZ);
			const MySparseMx &ZY = GetSparsePost(PairIndexZY);

			RelaxFlat_XZ_ZY(XZ, ZY, wZ, Post);
#if 0//TRACE
			LogFlatMx("ConsPair after RelaxFlat_XZ_ZY", Post, LX, LY);
#endif
			}
		else if (SeqIndexZ > SeqIndexX && SeqIndexZ > SeqIndexY)
			{
			uint PairIndexXZ = GetPairIndex(SeqIndexX, SeqIndexZ);
			uint PairIndexYZ = GetPairIndex(SeqIndexY, SeqIndexZ);

			const MySparseMx &XZ = GetSparsePost(PairIndexXZ);
			const MySparseMx &YZ = GetSparsePost(PairIndexYZ);

			RelaxFlat_XZ_YZ(XZ, YZ, wZ, Post);
#if 0//TRACE
			LogFlatMx("ConsPair after RelaxFlat_XZ_YZ", Post, LX, LY);
#endif
			}
		else
			asserta(false);
		}

	MySparseMx &UpdatedSparsePostXY = GetUpdatedSparsePost(PairIndex);

#if 0//TRACE
	LogFlatMx("Final post before update", Post, LX, LY);
#endif
	UpdatedSparsePostXY.UpdateFromPost(SparsePostXY, Post, SeqCount);
	UpdatedSparsePostXY.m_X = SparsePostXY.m_X;
	UpdatedSparsePostXY.m_Y = SparsePostXY.m_Y;
	myfree(Post);

#if 0//TRACE
	Log("\nBefore:");
	SparsePostXY.LogMe();
	Log("\nUpdated:");
	UpdatedSparsePostXY.LogMe();
#endif
	}
