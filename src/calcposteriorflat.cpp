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

	uint LX = GetL(SeqIndexX);
	uint LY = GetL(SeqIndexY);
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

	float *Fwd = AllocFB(LX, LY);
	float *Bwd = AllocFB(LX, LY);

	//CalcFwdFlat(X, LX, Y, LY, Fwd);
	//CalcBwdFlat(X, LX, Y, LY, Bwd);

	uint GSIX = GetGSIByLabel(LabelX);
	uint GSIY = GetGSIByLabel(LabelY);
	uint LX2 = GetSeqLengthByGSI(GSIX);
	uint LY2 = GetSeqLengthByGSI(GSIY);
	asserta(LX2 == LX);
	asserta(LY2 == LY);

	CalcFwdFlat_MPCFlat(GSIX, LX, GSIY, LY, Fwd);
	CalcBwdFlat_MPCFlat(GSIX, LX, GSIY, LY, Bwd);

	float *Post = AllocPost(LX, LY);
	CalcPostFlat(Fwd, Bwd, LX, LY, Post);
#if 0//TRACE
	LogFlatMxs("FwdFlat", Fwd, LX, LY);
	LogFlatMxs("BwdFlat", Bwd, LX, LY);
	LogFlatMx("PostFlat", Post, LX, LY);
#endif
	myfree(Fwd);
	myfree(Bwd);

#if 0//TRACE
	LogFlatMx1("Fwd", Fwd, LX, LY);
	LogFlatMx1("Bwd", Bwd, LX, LY);
	LogFlatMx("Post", Post, LX, LY);
#endif

	MySparseMx &SparsePost = GetSparsePost(PairIndex);
	SparsePost.FromPost(Post, LX, LY);

	const byte *X = GetBytePtr(SeqIndexX);
	const byte *Y = GetBytePtr(SeqIndexY);
	SparsePost.m_X = X;
	SparsePost.m_Y = Y;

#if 0//TRACE
	SparsePost.LogMe();
#endif

	float *DPRows = AllocDPRows(LX, LY);
	float Score = CalcAlnScoreFlat(Post, LX, LY, DPRows);
	myfree(Post);
	myfree(DPRows);

#if 0//TRACE
	string Path;
	char *TB = myalloc(char, (LX+1)*(LY+1));
	float Score2 = CalcAlnFlat(Post, LX, LY, DPRows, TB, Path);
	Log("Score=%.3g Score2=%.3g\n", Score, Score2);
	myfree(TB);
#endif

	float EA = Score/min(LX, LY);
#if 0//TRACE
	const char *LabelX = GetLabel(SeqIndexX);
	const char *LabelY = GetLabel(SeqIndexY);
	Log("Flat EA(%s, %s) = %.3g\n", LabelX, LabelY, EA);
#endif
	m_DistMx[SeqIndexX][SeqIndexY] = EA;
	m_DistMx[SeqIndexY][SeqIndexX] = EA;
	}
