#include "muscle.h"
#include "mpcflat.h"

void CalcPostFlat(const float *FlatFwd, const float *FlatBwd,
  uint LX, uint LY, float *Post)
	{
	float Total = CalcTotalProbFlat(FlatFwd, FlatBwd, LX, LY);

	uint IxFB = HMMSTATE_COUNT*((LY + 1) + 1); // M[1,1]
	uint IxPost = 0;
	for (uint i = 0; i < LX; ++i)
		{
		for (uint j = 0; j < LY; ++j)
			{
			float Score = FlatFwd[IxFB] + FlatBwd[IxFB] - Total;
			if (Score < MIN_SPARSE_SCORE)
				Post[IxPost++] = 0;
			else
				{
				float P = (Score >= LOG_ONE ? 1.0f : expf(Score));
				Post[IxPost++] = P;
				}
			IxFB += HMMSTATE_COUNT;
			}
		IxFB += HMMSTATE_COUNT;
		}
	}

void MPCFlat::CalcFwdFlat_MPCFlat(uint GSIX, uint LX,
  uint GSIY, uint LY, float *Flat)
	{
	const byte *X = GetByteSeqByGSI(GSIX);
	const byte *Y = GetByteSeqByGSI(GSIY);
	CalcFwdFlat(X, LX, Y, LY, Flat);
	}

void MPCFlat::CalcBwdFlat_MPCFlat(uint GSIX, uint LX,
  uint GSIY, uint LY, float *Flat)
	{
	const byte *X = GetByteSeqByGSI(GSIX);
	const byte *Y = GetByteSeqByGSI(GSIY);
	CalcBwdFlat(X, LX, Y, LY, Flat);
	}

void MPCFlat::CalcPosterior(uint PairIndex)
	{
	const pair<uint, uint> &Pair = GetPair(PairIndex);

	const uint SeqIndexX = Pair.first;
	const uint SeqIndexY = Pair.second;

	uint LX = m_MyInputSeqs->GetSeqLength(SeqIndexX);
	uint LY = m_MyInputSeqs->GetSeqLength(SeqIndexY);
	if (double(LX)*double(LY)*5 + 100 > double(INT_MAX))
		{
		ProgressLog("\nSequence length %u >%s\n",
		  LX, GetLabel(SeqIndexX));
		ProgressLog("Sequence length %u >%s\n",
		  LY, GetLabel(SeqIndexY));
		Die("HMM overflow, sequence lengths %u, %u (max ~21k)", LX, LY);
		}

	const string LabelX = string(m_MyInputSeqs->GetLabel(SeqIndexX));
	const string LabelY = string(m_MyInputSeqs->GetLabel(SeqIndexY));

	uint GSIX = GetGSIByLabel(LabelX);
	uint GSIY = GetGSIByLabel(LabelY);

	uint LX2 = GetSeqLengthByGSI(GSIX);
	uint LY2 = GetSeqLengthByGSI(GSIY);
	asserta(LX2 == LX);
	asserta(LY2 == LY);

	float *Post = CalcPost(LabelX, LabelY);

	MySparseMx &SparsePost = GetSparsePost(PairIndex);
	SparsePost.FromPost(Post, LX, LY);

	const byte *X = GetBytePtr(SeqIndexX);
	const byte *Y = GetBytePtr(SeqIndexY);
	SparsePost.m_X = X;
	SparsePost.m_Y = Y;

	float *DPRows = AllocDPRows(LX, LY);
	float Score = CalcAlnScoreFlat(Post, LX, LY, DPRows);
	myfree(Post);
	myfree(DPRows);

	float EA = Score/min(LX, LY);
	m_DistMx[SeqIndexX][SeqIndexY] = EA;
	m_DistMx[SeqIndexY][SeqIndexX] = EA;
	}
