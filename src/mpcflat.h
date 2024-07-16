#pragma once

#include "multisequence.h"
#include "upgma5.h"
#include "tree.h"
#include "treeperm.h"
#include "mysparsemx.h"
#include "derep.h"
#include "clustalweights.h"
#include <unordered_map>

static const uint DEFAULT_CONSISTENCY_ITERS_FLAT = 2;
static const uint DEFAULT_REFINE_ITERS_FLAT = 100;

// Multi-threaded ProbCons, flat memory layout
class MPCFlat
	{
public:
	MultiSequence *m_OriginalInputSeqs = 0;

// m_MyInputSeqs are unique seqs after derep of m_OriginalInputSeqs
	MultiSequence *m_MyInputSeqs = 0;

	MultiSequence *m_MSA = 0;

	uint m_ConsistencyIterCount = DEFAULT_CONSISTENCY_ITERS_FLAT;
	uint m_RefineIterCount = DEFAULT_REFINE_ITERS_FLAT;
	TREEPERM m_TreePerm = TP_None;
	vector<string> m_Labels;
	unordered_map<string, uint> m_LabelToIndex;
	UPGMA5 m_Upgma5;
	Tree m_GuideTree;
	vector<MultiSequence *> m_ProgMSAs;
	Derep m_D;

	ClustalWeights m_CW;
	vector<float> m_Weights;

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

public:
	virtual void CalcFwdFlat_MPCFlat(uint GSIX, uint LX,
	  uint GSIY, uint LY, float *Flat);
	virtual void CalcBwdFlat_MPCFlat(uint GSIX, uint LX,
	  uint GSIY, uint LY, float *Flat);

public:
	void AllocPairCount(uint SeqCount);
	void FreeProgMSAs();
	void FreeSparsePosts();
	void InitSeqs(MultiSequence *InputSeqs);
	void InitPairs();
	void InitDistMx();
	void CalcPosterior(uint PairIndex);
	void CalcPosteriors();
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
	uint GetMyInputSeqIndex(const string &Label) const;
	const byte *GetBytePtr(uint SeqIndex) const;
	uint GetPairIndex(uint SMI1, uint SMI2) const;
	MySparseMx &GetSparsePost(uint PairIndex);
	MySparseMx &GetUpdatedSparsePost(uint PairIndex);
	void BuildPost(const MultiSequence &MSA1, const MultiSequence &MSA2,
	  float *Post);
	uint GetSeqLength(uint SeqIndex) const;
	const Sequence *GetSequence(uint SeqIndex) const;
	MultiSequence *AlignAlns(const MultiSequence &MSA1,
	  const MultiSequence &MSA2, float *ptrScore = 0);
	void GetLabelToMSASeqIndex(unordered_map<string, uint> &LabelToMSASeqIndex) const;
	void SortMSA_ByInputOrder();
	void SortMSA_ByGuideTree();
	void SortMSA();
	void InsertDupes(const unordered_map<string, vector<string> > &RepSeqLabelToDupeLabels);
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

void RelaxFlat_XZ_ZY(const MySparseMx &XZ, const MySparseMx &ZY,
  float WeightZ, float *Post);
void RelaxFlat_ZX_ZY(const MySparseMx &XZ, const MySparseMx &YZ,
  float WeightZ, float *Post);
void RelaxFlat_XZ_YZ(const MySparseMx &XZ, const MySparseMx &YZ,
  float WeightZ, float *Post);
