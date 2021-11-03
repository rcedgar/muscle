#include "muscle.h"

void GetAllPairs(uint Count, 
  vector<uint> &Indexes1, vector<uint> &Indexes2)
	{
	Indexes1.clear();
	Indexes2.clear();
	for (uint i = 0; i < Count; ++i)
		{
		for (uint j = i + 1; j < Count; ++j)
			{
			Indexes1.push_back(i);
			Indexes2.push_back(j);
			}
		}
	}

void GetAllPairs(uint Count1, uint Count2,
  vector<uint> &Indexes1, vector<uint> &Indexes2)
	{
	Indexes1.clear();
	Indexes2.clear();
	for (uint i = 0; i < Count1; ++i)
		{
		for (uint j = 0; j < Count2; ++j)
			{
			Indexes1.push_back(i);
			Indexes2.push_back(j);
			}
		}
	}

void GetPairs(uint Count1, uint Count2, uint TargetPairCount,
  vector<uint> &Indexes1, vector<uint> &Indexes2)
	{
	Indexes1.clear();
	Indexes2.clear();

	uint AllPairCount = Count1*Count2;
	if (TargetPairCount == UINT_MAX || AllPairCount < TargetPairCount*3/2)
		{
		GetAllPairs(Count1, Count2, Indexes1, Indexes2);
		return;
		}

	set<pair<uint, uint> > PairSet;
	const uint MaxCounter = TargetPairCount*10;
	uint Counter = 0;
	while (Counter++ < MaxCounter && (uint) SIZE(PairSet) < TargetPairCount)
		{
		uint i = randu32()%Count1;
		uint j = randu32()%Count2;
		if (i == j)
			continue;
		pair<uint, uint> Pair(i, j);
		PairSet.insert(Pair);
		}

	uint PairCount = SIZE(PairSet);
	asserta(PairCount > TargetPairCount/2);
	for (set<pair<uint, uint> >::const_iterator p = PairSet.begin();
		p != PairSet.end(); ++p)
		{
		uint Index1 = p->first;
		uint Index2 = p->second;
		Indexes1.push_back(Index1);
		Indexes2.push_back(Index2);
		}
	}
