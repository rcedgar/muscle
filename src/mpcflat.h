#pragma once

#include "multisequence.h"
#include "upgma5.h"
#include "tree.h"
#include "treeperm.h"
#include "mysparsemx.h"

// Multi-threaded ProbCons
class MPCFlat
	{
public:
	MultiSequence *m_InputSeqs = 0;
	MultiSequence *m_MSA = 0;

	uint m_ConsistencyIterCount = DEFAULT_CONSISTENCY_ITERS;
	uint m_RefineIterCount = DEFAULT_REFINE_ITERS;
	TREEPERM m_TreePerm = TP_None;
	vector<string> m_Labels;
	map<string, uint> m_LabelToIndex;
	UPGMA5 m_Upgma5;
	Tree m_GuideTree;
	vector<MultiSequence *> m_ProgMSAs;

	vector<vector<float> > m_DistMx;
	vector<pair<uint, uint> > m_Pairs;
	map<pair<uint, uint>, uint> m_PairToIndex;
	vector<uint> m_JoinIndexes1;
	vector<uint> m_JoinIndexes2;

// Per-pair
	vector<MySparseMx *> m_SparsePosts1;
	vector<MySparseMx *> m_SparsePosts2;
	vector<MySparseMx *> *m_ptrSparsePosts = &m_SparsePosts1;
	vector<MySparseMx *> *m_ptrUpdatedSparsePosts = &m_SparsePosts2;

public:
	~MPCFlat()
		{
		Clear();
		}

	void Clear();
	void Run(MultiSequence *InputSeqs);
	uint GetSeqCount() const;
	void Run_Super4(MultiSequence *InputSeqs);

private:
	void AllocPairCount(uint SeqCount);
	void FreeProgMSAs();
	void FreeSparsePosts();
	uint GetL(uint SeqIndex) const;
	void InitSeqs(MultiSequence *InputSeqs);
	void InitPairs();
	void InitDistMx();
	void CalcPosteriors();
	void CalcPosterior(uint PairIndex);
	void Consistency();
	void ConsIter(uint Iter);
	void ConsPair(uint PairIndex);
	void CalcGuideTree();
	void CalcGuideTree_RandomChain();
	void CalcJoinOrder();
	void ProgressiveAlign();
	void Refine();
	void RefineIter();
	void ProgAln(uint JoinIndex);
	const pair<uint, uint> &GetPair(uint PairIndex) const;
	const char *GetLabel(uint SeqIndex) const;
	const byte *GetBytePtr(uint SeqIndex) const;
	uint GetPairIndex(uint SMI1, uint SMI2) const;
	MySparseMx &GetSparsePost(uint PairIndex);
	MySparseMx &GetUpdatedSparsePost(uint PairIndex);
	MultiSequence *AlignAlns(const MultiSequence &MSA1, const MultiSequence &MSA2);
	void BuildPost(const MultiSequence &MSA1, const MultiSequence &MSA2,
	  float *Post);
	uint GetSeqLength(uint SeqIndex) const;
	const Sequence *GetSequence(uint SeqIndex) const;
	};

float *AllocFB(uint LX, uint LY);
float *AllocPost(uint LX, uint LY);
float *AllocDPRows(uint LX, uint LY);
char *AllocTB(uint LX, uint LY);

float CalcTotalProbFlat(const float *FlatFwd, const float *FlatBwd, uint LX, uint LY);
void CalcFwdFlat(const byte *X, uint LX, const byte *Y, uint LY, float *Flat);
void CalcBwdFlat(const byte *X, uint LX, const byte *Y, uint LY, float *Flat);
float CalcAlnScoreFlat(const float *Post, uint LX, uint LY, float *DPRows);
float CalcAlnFlat(const float *Post, uint LX, uint LY,
  float *DPRows, char *TB, string &Path);

void RelaxFlat_XZ_ZY(const MySparseMx &XZ, const MySparseMx &ZY, float *Post);
void RelaxFlat_ZX_ZY(const MySparseMx &XZ, const MySparseMx &YZ, float *Post);
void RelaxFlat_XZ_YZ(const MySparseMx &XZ, const MySparseMx &YZ, float *Post);
