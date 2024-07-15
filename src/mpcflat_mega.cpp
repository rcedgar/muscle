// Created by Igor Tolstoy 2024-06-19.
// Updated by Igor Tolstoy & Robert Edgar

#include "muscle.h"
#include "mpcflat_mega.h"

//// SMI = Sparse Matrix (posterior) Index
////   same as index into m_InputSeqs << NOT WHEN PROFILE-PROFILE SUBSAMPLED
//const vector<vector<byte> > &MPCFlat_mega::GetProfile(uint SMI) const
//	{
//	asserta(SMI < SIZE(*m_ProfilePtrVec));
//	return *(*m_ProfilePtrVec)[SMI];
//	}
//
//void MPCFlat_mega::CalcPosterior(uint PairIndex)
//    {
//    const pair<uint, uint> &Pair = GetPair(PairIndex);
//
//    const uint SeqIndexX = Pair.first;
//    const uint SeqIndexY = Pair.second;
//
//    uint LX = GetL(SeqIndexX);
//    uint LY = GetL(SeqIndexY);
//    if (double(LX)*double(LY)*5 + 100 > double(INT_MAX))
//        {
//        ProgressLog("\nSequence length %u >%s\n",
//          LX, GetLabel(SeqIndexX));
//        ProgressLog("Sequence length %u >%s\n",
//          LY, GetLabel(SeqIndexY));
//        Die("HMM overflow, sequence lengths %u, %u (max ~21k)", LX, LY);
//        }
//
//	const vector<vector<byte> > &ProfileX = GetProfile(SeqIndexX);
//	const vector<vector<byte> > &ProfileY = GetProfile(SeqIndexY);
//
//	const uint PLX = SIZE(ProfileX);
//	const uint PLY = SIZE(ProfileY);
//
//	asserta(PLX == LX);
//	asserta(PLY == LY);
//
//    float *Fwd = AllocFB(LX, LY);
//    float *Bwd = AllocFB(LX, LY);
//
//    const byte *X = GetBytePtr(SeqIndexX);
//    const byte *Y = GetBytePtr(SeqIndexY);
//
////    CalcFwdFlat(X, LX, Y, LY, Fwd);
////    CalcBwdFlat(X, LX, Y, LY, Bwd);
//
//    Mega::CalcFwdFlat_mega(ProfileX, ProfileY, Fwd);
//    Mega::CalcBwdFlat_mega(ProfileX, ProfileY, Bwd);
//    
//    float *Post = AllocPost(LX, LY);
//    CalcPostFlat(Fwd, Bwd, LX, LY, Post);
//#if 0//TRACE
//    LogFlatMxs("FwdFlat", Fwd, LX, LY);
//    LogFlatMxs("BwdFlat", Bwd, LX, LY);
//    LogFlatMx("PostFlat", Post, LX, LY);
//#endif
//    myfree(Fwd);
//    myfree(Bwd);
//
//#if 0//TRACE
//    LogFlatMx1("Fwd", Fwd, LX, LY);
//    LogFlatMx1("Bwd", Bwd, LX, LY);
//    LogFlatMx("Post", Post, LX, LY);
//#endif
//
//    MySparseMx &SparsePost = GetSparsePost(PairIndex);
//    SparsePost.FromPost(Post, LX, LY);
//    SparsePost.m_X = X;
//    SparsePost.m_Y = Y;
//
//#if 0//TRACE
//    SparsePost.LogMe();
//#endif
//
//    float *DPRows = AllocDPRows(LX, LY);
//    float Score = CalcAlnScoreFlat(Post, LX, LY, DPRows);
//    myfree(Post);
//    myfree(DPRows);
//
//#if 0//TRACE
//    string Path;
//    char *TB = myalloc(char, (LX+1)*(LY+1));
//    float Score2 = CalcAlnFlat(Post, LX, LY, DPRows, TB, Path);
//    Log("Score=%.3g Score2=%.3g\n", Score, Score2);
//    myfree(TB);
//#endif
//
//    float EA = Score/min(LX, LY);
//#if 0//TRACE
//    const char *LabelX = GetLabel(SeqIndexX);
//    const char *LabelY = GetLabel(SeqIndexY);
//    Log("Flat EA(%s, %s) = %.3g\n", LabelX, LabelY, EA);
//#endif
//    m_DistMx[SeqIndexX][SeqIndexY] = EA;
//    m_DistMx[SeqIndexY][SeqIndexX] = EA;
//    }
//
//void MPCFlat_mega::Run(MultiSequence *InputSeqs,
//  const vector<const vector<vector<byte> > *> &ProfilePtrVec)
//	{
//	m_ProfilePtrVec = &ProfilePtrVec;
//	m_D.m_Disable = true;
//	MPCFlat::Run(InputSeqs);
//	}

void MPCFlat_mega::CalcFwdFlat_MPCFlat(uint GSIX, uint LX,
  uint GSIY, uint LY, float *Flat)
	{
	const vector<vector<byte> > &ProfileX = *Mega::GetProfileByGSI(GSIX);
	const vector<vector<byte> > &ProfileY = *Mega::GetProfileByGSI(GSIY);
	asserta(SIZE(ProfileX) == LX);
	asserta(SIZE(ProfileY) == LY);
    Mega::CalcFwdFlat_mega(ProfileX, ProfileY, Flat);
	}

void MPCFlat_mega::CalcBwdFlat_MPCFlat(uint GSIX, uint LX,
  uint GSIY, uint LY, float *Flat)
	{
	const vector<vector<byte> > &ProfileX = *Mega::GetProfileByGSI(GSIX);
	const vector<vector<byte> > &ProfileY = *Mega::GetProfileByGSI(GSIY);
	asserta(SIZE(ProfileX) == LX);
	asserta(SIZE(ProfileY) == LY);
    Mega::CalcBwdFlat_mega(ProfileX, ProfileY, Flat);
	}
