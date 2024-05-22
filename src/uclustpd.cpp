#include "muscle.h"
#include "uclustpd.h"
#include "objmgr.h"
#include "pathinfo.h"
#include "omplock.h"
#include "sort.h"

float ViterbiFastMem(XDPMem &Mem, const byte *A, uint LA,
  const byte *B, uint LB, PathInfo &PI);
void LogAln(const byte *X, uint LX, const byte *Y, uint LY, const PathInfo &PI);
void MakeAlnRows(const byte *XSeq, uint LX,
  const byte *YSeq, uint LY, const PathInfo &PI,
  string &RowX, string &RowY);
double GetProtDist(const char *Q, const char *T, uint ColCount);

static vector<XDPMem *> g_Mems;
XDPMem &GetDPMem()
	{
	LOCK();
	if (SIZE(g_Mems) == 0)
		{
		uint n = GetRequestedThreadCount();
		asserta(n > 0);
		for (uint i = 0; i < n; ++i)
			g_Mems.push_back(new XDPMem);
		}
	uint ThreadIndex = GetThreadIndex();
	asserta(ThreadIndex < SIZE(g_Mems));
	UNLOCK();
	return *g_Mems[ThreadIndex];
	}

const char *UClustPD::GetLabel(uint SeqIndex) const
	{
	const char *Label = m_InputSeqs->GetLabel(SeqIndex);
	return Label;
	}

const byte *UClustPD::GetByteSeq(uint SeqIndex, uint &L) const
	{
	const byte *Seq = m_InputSeqs->GetByteSeq(SeqIndex, L);
	return Seq;
	}

double UClustPD::GetProtDistPair(uint SeqIndexi, uint SeqIndexj, string *Path)
	{
	uint Li;
	const byte *Seqi = GetByteSeq(SeqIndexi, Li);
	const char *Labeli = GetLabel(SeqIndexi);

	uint Lj;
	const byte *Seqj = GetByteSeq(SeqIndexj, Lj);
	const char *Labelj = GetLabel(SeqIndexj);

	PathInfo *PI = ObjMgr::GetPathInfo();
	XDPMem &Mem = GetDPMem();
	ViterbiFastMem(Mem, Seqi, Li, Seqj, Lj, *PI);

	string RowX, RowY;
	MakeAlnRows(Seqi, Li, Seqj, Lj, *PI, RowX, RowY);

	const uint ColCount = PI->GetColCount();
	asserta(SIZE(RowX) == ColCount);
	asserta(SIZE(RowY) == ColCount);
	double dij = ::GetProtDist(RowX.c_str(), RowY.c_str(), ColCount);
	if (Path != 0)
		PI->GetPathStr(*Path);
	ObjMgr::Down(PI);
	return dij;
	}

uint UClustPD::SearchAll(uint SeqIndex)
	{
	const int N = (int) GetSubsetSize();
	uint HitCount = 0;
#pragma omp parallel for num_threads(m_ThreadCount)
	for (int i = 0; i < N; ++i)
		{
		uint SeqIndex2 = (*m_SubsetSeqIndexes)[i];
		if (SeqIndex2 == SeqIndex)
			continue;
		double d = GetProtDistPair(SeqIndex, SeqIndex2);
		if (d < m_MaxPD)
			{
			LOCK();
			++HitCount;
			UNLOCK();
			}
		}
	return HitCount;
	}

uint UClustPD::Search(uint SeqIndex, const vector<uint> &Centroids,
  double &BestDist)
	{
	const int N = (int) SIZE(Centroids);
	uint BestCentroidIndex = UINT_MAX;
	BestDist = FLT_MAX;
#pragma omp parallel for num_threads(m_ThreadCount)
	for (int i = 0; i < N; ++i)
		{
		uint CentroidIndex = Centroids[i];
		asserta(CentroidIndex < SIZE(m_CentroidSeqIndexes));
		uint CentroidSeqIndex = m_CentroidSeqIndexes[CentroidIndex];
		//string Path;
		double d = GetProtDistPair(SeqIndex, CentroidSeqIndex);
		if (d > m_MaxPD)
			continue;
		LOCK();
		if (d < BestDist)
			{
			BestDist = d;
			BestCentroidIndex = CentroidIndex;
			}
		UNLOCK();
		}
	return BestCentroidIndex;
	}

void UClustPD::CentroidsToFasta(const string &FileName) const
	{
	if (FileName == "")
		return;
	FILE *f = CreateStdioFile(FileName);
	const uint ClusterCount = GetClusterCount();
	asserta(SIZE(m_CentroidSeqIndexes) == ClusterCount);
	for (uint ClusterIndex = 0; ClusterIndex < ClusterCount;
	  ++ClusterIndex)
		{
		uint CentroidSeqIndex = m_CentroidSeqIndexes[ClusterIndex];
		const Sequence *Seq = m_InputSeqs->GetSequence(CentroidSeqIndex);
		Seq->WriteMFA(f);
		}
	CloseStdioFile(f);
	}

void UClustPD::Run(MultiSequence &InputSeqs, 
  const vector<uint> &SubsetSeqIndexes, double MaxPD)
	{
	Clear();
	m_ThreadCount = GetRequestedThreadCount();
	m_MaxPD = MaxPD;
	m_InputSeqs = &InputSeqs;
	m_SubsetSeqIndexes = &SubsetSeqIndexes;

	const uint SubsetSize = GetSubsetSize();
	for (uint i = 0; i < SubsetSize; ++i)
		m_PendingSubsetIndexes.push_back(i);

	m_SubsetIndexToCentroidIndex.resize(SubsetSize, UINT_MAX);
	m_SubsetIndexToDist.resize(SubsetSize, 0);
	uint Iter = 0;
	for (;;)
		{
		++Iter;
		uint PendingCount = SIZE(m_PendingSubsetIndexes);
		uint DoneCount = SubsetSize - PendingCount;
		double DonePct = GetPct(DoneCount, SubsetSize);
		Progress("Iter [%u] %u clusters, (assigned %.3g%%, remaining %s)\n",
		  Iter, SIZE(m_CentroidSeqIndexes), DonePct, IntToStr(PendingCount));
		uint BeginSize = SIZE(m_PendingSubsetIndexes);
		if (BeginSize == 0)
			break;
		vector<uint> NewSeeds;
		vector<list<uint>::iterator> DoneVec1;
		uint PendingIndex = 0;
		for (list<uint>::iterator p = m_PendingSubsetIndexes.begin();
		  p != m_PendingSubsetIndexes.end(); ++p)
			{
			ProgressStep(PendingIndex++, PendingCount,
			  "New seeds (%u)", SIZE(NewSeeds));
			uint SubsetIndex = *p;
			asserta(SubsetIndex < SIZE(*m_SubsetSeqIndexes));
			uint SeqIndex = (*m_SubsetSeqIndexes)[SubsetIndex];

			double d;
			uint CentroidIndex = Search(SeqIndex, NewSeeds, d);
			if (CentroidIndex == UINT_MAX)
				{
				uint CentroidIndex = SIZE(m_CentroidSeqIndexes);
				NewSeeds.push_back(CentroidIndex);
				m_SubsetIndexToCentroidIndex[SubsetIndex] = CentroidIndex;
				m_SubsetIndexToDist[SubsetIndex] = 0;
				m_CentroidSeqIndexes.push_back(SeqIndex);

				vector<uint> NewVec;
				NewVec.push_back(SubsetIndex);
				m_CentroidIndexToMemberSubsetIndexes.push_back(NewVec);
				asserta(SIZE(m_CentroidSeqIndexes) ==
				  SIZE(m_CentroidIndexToMemberSubsetIndexes));
				DoneVec1.push_back(p);
				}
			if (SIZE(NewSeeds) >= m_ThreadCount)
				{
				ProgressStep(PendingCount-1, PendingCount, "New seeds (%u)",
				  SIZE(NewSeeds));
				break;
				}
			}
		asserta(!DoneVec1.empty());
		if (DoneVec1.size() == m_PendingSubsetIndexes.size())
			{
			m_PendingSubsetIndexes.clear();
			break;
			}
		for (uint i = 0; i < SIZE(DoneVec1); ++i)
			m_PendingSubsetIndexes.erase(DoneVec1[i]);

		PendingCount = SIZE(m_PendingSubsetIndexes);
		if (PendingCount == 0)
			break;

		PendingIndex = 0;
		uint NewHitCount = 0;
		vector<list<uint>::iterator> DoneVec2;
		for (list<uint>::iterator p = m_PendingSubsetIndexes.begin();
		  p != m_PendingSubsetIndexes.end(); ++p)
			{
			ProgressStep(PendingIndex++, PendingCount, "Align (%u)", NewHitCount);
			uint SubsetIndex = *p;
			asserta(SubsetIndex < SIZE(*m_SubsetSeqIndexes));
			uint SeqIndex = (*m_SubsetSeqIndexes)[SubsetIndex];

			double d;
			uint CentroidIndex = Search(SeqIndex, NewSeeds, d);
			if (CentroidIndex != UINT_MAX)
				{
				asserta(CentroidIndex < SIZE(m_CentroidSeqIndexes));
				uint CentroidSeqIndex = m_CentroidSeqIndexes[CentroidIndex];
				const char *LabelQ = GetLabel(SeqIndex);
				const char *LabelC = GetLabel(CentroidSeqIndex);
				m_SubsetIndexToCentroidIndex[SubsetIndex] = CentroidIndex;
				m_SubsetIndexToDist[SubsetIndex] = d;
				asserta(CentroidIndex < SIZE(m_CentroidIndexToMemberSubsetIndexes));
				m_CentroidIndexToMemberSubsetIndexes[CentroidIndex].push_back(SubsetIndex);
				++NewHitCount;
				DoneVec2.push_back(p);
				}
			}
		
		if (DoneVec2.size() == m_PendingSubsetIndexes.size())
			{
			m_PendingSubsetIndexes.clear();
			break;
			}
		for (uint i = 0; i < SIZE(DoneVec2); ++i)
			m_PendingSubsetIndexes.erase(DoneVec2[i]);
		uint EndSize = SIZE(m_PendingSubsetIndexes);
		asserta(EndSize < BeginSize);
		}
	}

uint UClustPD::GetClusterSize(uint ClusterIndex) const
	{
	asserta(ClusterIndex < SIZE(m_CentroidIndexToMemberSubsetIndexes));
	uint Size = SIZE(m_CentroidIndexToMemberSubsetIndexes[ClusterIndex]);
	return Size;
	}

void UClustPD::ToTsv(const string &FileName) const
	{
	if (FileName == "")
		return;
	FILE *f = CreateStdioFile(FileName);
	ToTsv(f);
	CloseStdioFile(f);
	}

void UClustPD::ToTsv(FILE *f) const
	{
	if (f == 0)
		return;

	const uint CentroidCount = SIZE(m_CentroidSeqIndexes);
	asserta(SIZE(m_CentroidIndexToMemberSubsetIndexes) == CentroidCount);
	const uint SubsetSize = GetSubsetSize();
	asserta(SIZE(m_SubsetIndexToCentroidIndex) == SubsetSize);

	uint DoneCount = 0;
	for (uint CentroidIndex = 0; CentroidIndex < CentroidCount;
	  ++CentroidIndex)
		{
		uint CentroidSeqIndex = m_CentroidSeqIndexes[CentroidIndex];
		const char *CentroidLabel = GetLabel(CentroidSeqIndex);
		const vector<uint> &MemberSubsetIndexes = 
		  m_CentroidIndexToMemberSubsetIndexes[CentroidIndex];
		const uint n = SIZE(MemberSubsetIndexes);
		asserta(n > 0);
		for (uint i = 0; i < n; ++i)
			{
			uint SubsetIndex = MemberSubsetIndexes[i];
			uint SeqIndex = (*m_SubsetSeqIndexes)[SubsetIndex];
			asserta(SubsetIndex < SIZE(m_SubsetIndexToDist));
			double d = m_SubsetIndexToDist[SubsetIndex];
			const char *Label = GetLabel(SeqIndex);
			fprintf(f, "%u\t%s\n", CentroidIndex, Label);
			}
		DoneCount += n;
		}
	asserta(DoneCount == SubsetSize);
	}

void UClustPD::GetClusterMFA(uint ClusterIndex, MultiSequence &MFA) const
	{
	MFA.Clear();
	const uint ClusterSize = GetClusterSize(ClusterIndex);
	asserta(ClusterIndex < SIZE(m_CentroidIndexToMemberSubsetIndexes));
	const vector<uint> &SubsetIndexes =
	  m_CentroidIndexToMemberSubsetIndexes[ClusterIndex];
	asserta(SIZE(SubsetIndexes) == ClusterSize);
	for (uint i = 0; i < ClusterSize; ++i)
		{
		uint SubsetIndex = SubsetIndexes[i];
		asserta(SubsetIndex < SIZE(*m_SubsetSeqIndexes));
		uint SeqIndex = (*m_SubsetSeqIndexes)[SubsetIndex];
		const Sequence *Seq = m_InputSeqs->GetSequence(SeqIndex);
		MFA.AddSequence(Seq, false);
		}
	}

void UClustPD::GetClusterMFAs(vector<MultiSequence *> &MFAs) const
	{
	MFAs.clear();
	const uint ClusterCount = GetClusterCount();
	for (uint ClusterIndex = 0; ClusterIndex < ClusterCount;
	  ++ClusterIndex)
		{
		MultiSequence *ClusterMFA = new MultiSequence;
		GetClusterMFA(ClusterIndex, *ClusterMFA);
		MFAs.push_back(ClusterMFA);
		}
	}

void UClustPD::GetClusterSizes(vector<uint> &Sizes) const
	{
	Sizes.clear();
	const uint ClusterCount = GetClusterCount();
	for (uint i = 0; i < ClusterCount; ++i)
		{
		uint Size = GetClusterSize(i);
		Sizes.push_back(Size);
		}
	}

void UClustPD::LogStats() const
	{
	uint ClusterCount = GetClusterCount();
	uint SubsetSize = GetSubsetSize();
	double AvgSize = SubsetSize / ClusterCount;
	vector<uint> ClusterSizes;
	GetClusterSizes(ClusterSizes);
	uint SingletonCount = 0;
	for (uint i = 0; i < ClusterCount; ++i)
		if (ClusterSizes[i] == 1)
			++SingletonCount;

	vector<uint> Order(ClusterCount);
	QuickSortOrderDesc(ClusterSizes.data(), ClusterCount, Order.data());
	uint MedianSize = ClusterSizes[Order[ClusterCount/2]];

	Log("\n");
	ProgressLog("%u seqs, %u clusters, avg size %.1f, median %u, singletons %u\n",
	  SubsetSize, ClusterCount, AvgSize, MedianSize, SingletonCount);

	for (uint i = 0; i < min(10u, ClusterCount); ++i)
		{
		uint k = Order[i];
		//         "   Avg. size  %.1f\n", AvgSize);
		ProgressLog("    Cluster  [%5u]   size %u\n", k, ClusterSizes[k]);
		}
	}

void cmd_uclustpd()
	{
	const string &InputFileName = opt(uclustpd);
	asserta(optset_maxpd);
	double MaxPD = opt(maxpd);
	if (optset_output)
		Die("Use -tsvout not -output");
	FILE *fOut = CreateStdioFile(opt(tsvout));

	MultiSequence Input;
	Input.FromFASTA(InputFileName, true);
	const uint SeqCount = Input.GetSeqCount();
	ProgressLog("%u seqs, maxpd %.2f\n", SeqCount, MaxPD);

	bool IsNucleo = Input.GuessIsNucleo();
	SetAlphab(IsNucleo);

	vector<uint> AllSeqIndexes;
	for (uint i = 0; i < SeqCount; ++i)
		AllSeqIndexes.push_back(i);

	UClustPD UD;
	UD.Run(Input, AllSeqIndexes, MaxPD);

	UD.ToTsv(fOut);
	CloseStdioFile(fOut);

	UD.LogStats();
	}
