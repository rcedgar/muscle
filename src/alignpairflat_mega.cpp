//
//  alignpairflat_mega.cpp
//  muscle5
//
//  Created by Igor Tolstoy on 6/17/24.
//

#include "alignpairflat_mega.hpp"

#include "muscle.h"
#include "mega.h"

float AlignPairFlat_mega_SparsePost(const Mega * M,
string &Path, MySparseMx *SparsePost, uint index_X, uint index_Y )
    {
    InitProbcons();

    
    asserta(SIZE(M->m_Seqs) > index_X);
    asserta(SIZE(M->m_Seqs) > index_Y);
    
uint L1 = SIZE(M->m_Seqs[index_X]);
    uint L2 = SIZE(M->m_Seqs[index_Y]);

//    uint L1 = Seq1->GetLength();
//    uint L2 = Seq2->GetLength();
    asserta(L1 > 0);
    asserta(L2 > 0);

//    const byte *ByteSeq1 = Seq1->GetBytePtr();
//    const byte *ByteSeq2 = Seq2->GetBytePtr();

    float *Fwd = AllocFB(L1, L2);
    float *Bwd = AllocFB(L1, L2);
    float *Post = AllocPost(L1, L2);

//    void CalcFwdFlat_mega(const Mega &M, uint ProfileIdxX, uint ProfileIdxY, float *Flat)
    CalcFwdFlat_mega(*M, index_X, index_Y, Fwd);

//  CalcFwdFlat(ByteSeq1, L1, ByteSeq2, L2, Fwd);
    
//    void CalcBwdFlat_mega(const Mega &M,uint ProfileIdxX, uint ProfileIdxY, float *Flat)
    CalcBwdFlat_mega(*M, index_X, index_Y, Bwd);
 
//    CalcBwdFlat(ByteSeq1, L1, ByteSeq2, L2, Bwd);
    CalcPostFlat(Fwd, Bwd, L1, L2, Post);
    delete Fwd;
    delete Bwd;

    float *DPRows = AllocDPRows(L1, L2);
    char *TB = AllocTB(L1, L2);
    float Score = CalcAlnFlat(Post, L1, L2, DPRows, TB, Path);
    if (SparsePost != 0)
        SparsePost->FromPost(Post, L1, L2);
    delete Post;
    delete DPRows;
    delete TB;

    asserta(L1 > 0 && L2 > 0);
    float EA = Score/min(L1, L2);
    return EA;
    }

float AlignPairFlat_mega(const Mega * M, string &Path, uint index_X,uint index_Y)
    {
    float EA = AlignPairFlat_mega_SparsePost(M, Path, 0, index_X, index_Y);
    return EA;
    }
