#pragma once

#include "multisequence.h"
#include "upgma5.h"
#include "tree.h"
#include "treeperm.h"

// Multi-threaded ProbCons
class MPC
	{
public:
	MultiSequence *m_InputSeqs = 0;
	MultiSequence *m_MSA = 0;

	uint m_ConsistencyIterCount = DEFAULT_CONSISTENCY_ITERS;
	uint m_RefineIterCount = DEFAULT_REFINE_ITERS;
	TREEPERM m_TreePerm = TP_None;
	vector<string> m_Labels;
	map<string, uint> m_LabelToIndex;
	UPGMA5 m_U5;
	Tree m_GuideTree;

	vector<vector<float> > m_DistMx;
	vector<vector<SparseMatrix *> > m_SparseMatrices;
	vector<pair<uint, uint> > m_Pairs;
	vector<uint> m_JoinIndexes1;
	vector<uint> m_JoinIndexes2;

	vector<vector<float> > m_AccAlnRows;
	vector<float> m_EAs;

public:
	void Run(MultiSequence *InputSeqs);
	uint GetSeqCount() const;
	void WriteAccAln(const string &FileName) const;

private:
	void InitSeqs(MultiSequence *InputSeqs);
	void InitPairs();
	void InitDistMx();
	void InitSparseMatrices();
	void CalcPosteriors();
	void CalcPosterior(uint SeqIndex1, uint SeqIndex2);
	void Consistency();
	void DoRelax(uint Iter);
	void CalcGuideTree();
	void CalcJoinOrder();
	void ProgressiveAlign();
	void Refine();
	float CalcAccAlnRow(uint PairIndex);
	void CalcAccAlnRows();
	void DeleteSparseMatrices();
	void WriteAccAln1(FILE *f, uint PairIndex) const;
	};
