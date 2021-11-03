#pragma once

#include "tree.h"

static const uint DEFAULT_TARGET_PAIR_COUNT = 2000;
static const uint DEFAULT_MAX_COARSE_SEQS = 500;

class PProg
	{
public:
	uint m_InputMSACount = 0;
	uint m_JoinCount = 0;
	uint m_NodeCount = 0;
	uint m_TargetPairCount = DEFAULT_TARGET_PAIR_COUNT;
	uint m_MaxCoarseSeqs = DEFAULT_MAX_COARSE_SEQS;

	map<string, uint> m_MSALabelToIndex;
	vector<string> m_MSALabels;
	vector<const MultiSequence *> m_MSAs;

	vector<uint> m_Pending;
	vector<vector<float> > m_ScoreMx;
	vector<vector<string> > m_PathMx;
	uint m_JoinIndex = 0;

	vector<uint> m_JoinMSAIndexes1;
	vector<uint> m_JoinMSAIndexes2;

public:
	void LoadMSAs(const vector<string> &FileNames, bool &IsNucleo);
	void SetMSAs(const vector<const MultiSequence *> &MSAs,
	  const vector<string> &MSALabels);

	void Run();
	void RunGuideTree(const Tree &GuideTree);
	void Run2(const vector<uint> &Indexes1,
	  const vector<uint> &Indexes2);
	void DeleteIndexesFromPending(uint Index1, uint Index2);
	void FindBestPair(uint &BestIndex1, uint &BestIndex2) const;
	void AlignAllInputPairs();
	void AlignAndJoin(uint Index1, uint Index2);
	void Join_ByPrecomputedPath(uint Index1, uint Index2);
	void AlignNewToPending();
	void LogPending(const string &s) const;
	const MultiSequence &GetMSA(uint Index) const;
	const string &GetMSALabel(uint Index) const;
	void SetMSA(uint Index, const MultiSequence &MSA);
	void SetMSALabel(uint Index, const string &Label);
	const MultiSequence &GetFinalMSA() const;
	void WriteGuideTree(const string &FileName) const;
	};

void MakeGuideTreeFromJoinOrder(const vector<uint> &Indexes1,
  const vector<uint> &Indexes2, const map<string, uint> &LabelToIndex,
  Tree &GuideTree);

void GetGuideTreeJoinOrder(const Tree &GuideTree,
  const map<string, uint> &LabelToIndex,
  vector<uint> &Indexes1, vector<uint> &Indexes2);

void ValidateJoinOrder(const vector<uint> &Indexes1,
  const vector<uint> &Indexes2);
