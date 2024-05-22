#include "muscle.h"
#include "uclustpd.h"
#include "sort.h"

void UClustPD::SelectCandidateGoodCentroids(vector<uint> &SubsetIndexes)
	{
	SubsetIndexes.clear();
	vector<uint> UnsortedSelectedIndexes;
	const uint ClusterCount = GetClusterCount();
	asserta(SIZE(m_CentroidIndexToMemberSubsetIndexes) == ClusterCount);
	vector<uint> SelectedHitCounts;
	for (uint ClusterIndex = 0; ClusterIndex < ClusterCount;
	  ++ClusterIndex)
		{
		ProgressStep(ClusterIndex, ClusterCount, "Select candidates");
		uint ClusterSize = GetClusterSize(ClusterIndex);
		if (ClusterSize < 8)
			continue;

		const vector<uint> &MemberSubsetIndexes =
		  m_CentroidIndexToMemberSubsetIndexes[ClusterIndex];
		uint n = uint(log2(ClusterSize)/2.0);
		asserta(n >= 1);
		set<uint> SampleSet;
		for (uint i = 0; i < 2*n; ++i)
			{
			uint r = randu32()%ClusterSize;
			uint SubsetIndex = MemberSubsetIndexes[r];
			SampleSet.insert(SubsetIndex);
			if (SIZE(SampleSet) == n)
				break;
			}

		uint TopHitCount = 0;
		uint TopSubsetIndex = UINT_MAX;
		for (set<uint>::const_iterator p = SampleSet.begin();
		  p != SampleSet.end(); ++p)
			{
			uint SubsetIndex = *p;
			asserta(SubsetIndex < SIZE(*m_SubsetSeqIndexes));
			uint SeqIndex = (*m_SubsetSeqIndexes)[SubsetIndex];
			uint HitCount = SearchAll(SubsetIndex);
			if (HitCount > TopHitCount)
				{
				TopHitCount = HitCount;
				TopSubsetIndex = SubsetIndex;
				}
			}
		Log("Cluster %u tophits %u size %u\n",
		  ClusterIndex, TopHitCount, ClusterSize);
		if (TopHitCount > ClusterSize)
			{
			UnsortedSelectedIndexes.push_back(TopSubsetIndex);
			SelectedHitCounts.push_back(TopHitCount);
			}
		else
			{
			uint CurrentCentroid =
			  m_CentroidIndexToMemberSubsetIndexes[ClusterIndex][0];
			UnsortedSelectedIndexes.push_back(CurrentCentroid);
			SelectedHitCounts.push_back(ClusterSize);
			}
		}
	const uint N = SIZE(UnsortedSelectedIndexes);
	uint *Order = myalloc(uint, N);
	QuickSortOrderDesc(SelectedHitCounts.data(), N, Order);
	for (uint k = 0; k < N; ++k)
		{
		uint i = Order[k];
		uint Index = UnsortedSelectedIndexes[i];
		uint HitCount = SelectedHitCounts[i];
		SubsetIndexes.push_back(Index);
		Log("Sorted %4u  hits %u\n", Index, HitCount);
		}
	myfree(Order);
	}

void cmd_uclustpd2()
	{
	const string &InputFileName = opt(uclustpd2);
	asserta(optset_maxpd);
	double MaxPD = opt(maxpd);
	if (optset_output || optset_tsvout)
		Die("Use -output1/2");

	MultiSequence Input;
	Input.FromFASTA(InputFileName, true);

	bool IsNucleo = Input.GuessIsNucleo();
	SetAlphab(IsNucleo);

	const uint SeqCount = Input.GetSeqCount();
	ProgressLog("%u sequences\n", SeqCount);

	vector<uint> AllSeqIndexes;
	for (uint i = 0; i < SeqCount; ++i)
		AllSeqIndexes.push_back(i);

	UClustPD UD;
	UD.Run(Input, AllSeqIndexes, MaxPD);

	UD.LogStats();
	UD.ToTsv(opt(output1));

	vector<uint> SIs;
	UD.SelectCandidateGoodCentroids(SIs);

	set<uint> SISet;
	vector<uint> AllSeqIndexes2;
	for (uint i = 0; i < SIZE(SIs); ++i)
		{
		uint Index = SIs[i];
		uint SeqIndex = AllSeqIndexes[Index];
		asserta(Index == SeqIndex);
		SISet.insert(Index);
		AllSeqIndexes2.push_back(SeqIndex);
		}

	for (uint i = 0; i < SeqCount; ++i)
		{
		if (SISet.find(i) == SISet.end())
			AllSeqIndexes2.push_back(i);
		}

	UD.Run(Input, AllSeqIndexes2, MaxPD);
	UD.LogStats();
	UD.ToTsv(opt(output2));
	UD.CentroidsToFasta(opt(centroids));
	}
