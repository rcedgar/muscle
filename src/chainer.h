#pragma once

#include "hspfinder.h"
#include <list>

const float BAD_SCORE = -9e9f;
struct HSPData;

// Bendpoint
struct BPData
	{
	uint Pos;
	bool IsLo;
	uint Index;

	void LogMe() const
		{
		Log("BP%s Pos %u Ix %u", (IsLo ? "lo" : "hi"), Pos, Index);
		}
	};

struct ChainData
	{
	uint LastHSPIndex;
	uint Ahi;
	uint Bhi;
	float Score;
	};

class Chainer
	{
public:
	const vector<const HSPData *> *m_HSPs;
	BPData *m_BPs;
	uint *m_PrevHSPIndexes;		// Predecessor in chain
	float *m_HSPIndexToChainScore;
	list<uint> m_Chains;		// Live HSP indexes

public:
	Chainer()
		{
		m_HSPs = 0;
		m_BPs = 0;
		m_PrevHSPIndexes = 0;
		m_HSPIndexToChainScore = 0;
		}

	~Chainer()
		{
		Clear();
		}

	void Clear()
		{
		m_HSPs = 0;
		m_Chains.clear();
		myfree(m_BPs);
		myfree(m_PrevHSPIndexes);
		myfree(m_HSPIndexToChainScore);
		m_BPs = 0;
		m_PrevHSPIndexes = 0;
		m_HSPIndexToChainScore = 0;
		}

	void Run(const vector<const HSPData *> &HSPs,
	  vector<const HSPData *> &Chain);
	void LogMe() const;
	void LogBPs() const;

public:
	static void LogHSPs(const vector<const HSPData *> &HSPs);
	static bool IsValidChain(const vector<const HSPData *> &HSPs);
	static void AssertValidChain(const vector<const HSPData *> &HSPs);
	static void LogChain(const vector<const HSPData *> &HSPs,
	  bool TestValid);
	static float GetChainScore(HSPData **HSPs, uint HSPCount);

private:
	void SetBPs();
	void SortBPs();
	uint FindBestChainLT(uint Ahi, uint Bhi);
	};
