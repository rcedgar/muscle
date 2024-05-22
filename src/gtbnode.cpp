#if 0

#include "muscle.h"
#include "gtbuilder.h"
#include "pathinfo.h"
#include "objmgr.h"
#include "omplock.h"

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
	if (SIZE(g_Mems) == 0)
		{
		uint n = GetRequestedThreadCount();
		asserta(n > 0);
		for (uint i = 0; i < n; ++i)
			g_Mems.push_back(new XDPMem);
		}
	uint ThreadIndex = GetThreadIndex();
	asserta(ThreadIndex < SIZE(g_Mems));
	return *g_Mems[ThreadIndex];
	}

const char *GTBNode::GetLabel(uint SeqIndex) const
	{
	const char *Label = m_Builder->m_Seqs->GetLabel(SeqIndex);
	return Label;
	}

const byte *GTBNode::GetByteSeq(uint SeqIndex, uint &L) const
	{
	const byte *Seq = m_Builder->m_Seqs->GetByteSeq(SeqIndex, L);
	return Seq;
	}

double GTBNode::GetProtDist(uint SeqIndexi, uint SeqIndexj)
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
	ObjMgr::Down(PI);
	return dij;
	}

void GTBNode::DoAll()
	{
	uint SeqCount = SIZE(m_SeqIndexes);
	asserta(SeqCount > 0);
	if (SeqCount < 3)
		return;
	vector<vector<float> > &DistMx = m_UPGMA.m_DistMx;
	DistMx.resize(SeqCount);
	m_UPGMA.m_Labels.clear();
	for (uint i = 0; i < SeqCount; ++i)
		{
		uint SeqIndex = m_SeqIndexes[i];
		const char *Label = m_Builder->m_Seqs->GetLabel(SeqIndex);
		m_UPGMA.m_Labels.push_back(Label);

		DistMx[i].resize(SeqCount, FLT_MAX);
		DistMx[i][i] = 0;
		}

	uint PairCount = (SeqCount*(SeqCount-1))/2;
	uint i = 1;
	uint j = 0;
	uint PairCounter = 0;
	uint ThreadCount = GetRequestedThreadCount();
#pragma omp parallel for num_threads(ThreadCount)
	for (int PairIndex = 0; PairIndex < (int) PairCount; ++PairIndex)
		{
		LOCK();
		ProgressStep(PairCounter++, PairCount, "All pairs");
		UNLOCK();
		uint SeqIndexi = m_SeqIndexes[i];
		uint SeqIndexj = m_SeqIndexes[j];
		double d = GetProtDist(SeqIndexi, SeqIndexj);
		DistMx[i][j] = float(d);
		DistMx[j][i] = float(d);

		++j;
		assert(j <= i);
		if (j == i)
			{
			++i;
			j = 0;
			}
		}
	m_UPGMA.Run(LINKAGE_Avg, m_Tree);
	}

void GTBNode::Run()
	{
	asserta(m_Builder != 0);
	asserta(!m_SeqIndexes.empty());
	const uint SeqCount = SIZE(m_SeqIndexes);
	if (SeqCount < m_Builder->m_MaxAll)
		{
		DoAll();
		return;
		}
	uint n = min(m_Builder->m_TargetSeedCount, SeqCount);
	set<uint> SeedSet;
	for (uint i = 0; i < n + 4; ++i)
		{
		uint r = randu32()%SeqCount;
		uint SeqIndex = m_SeqIndexes[i];
		SeedSet.insert(SeqIndex);
		if (SIZE(SeedSet) == n)
			break;
		}
	uint SeedCount = SIZE(SeedSet);
	for (set<uint>::const_iterator p = SeedSet.begin();
	  p != SeedSet.end(); ++p)
		m_SeedSeqIndexes.push_back(*p);
	asserta(SIZE(m_SeedSeqIndexes) == SeedCount);

	vector<vector<float> > &DistMx = m_UPGMA.m_DistMx;
	DistMx.resize(SeedCount);
	m_UPGMA.m_Labels.clear();
	for (uint i = 0; i < SeedCount; ++i)
		{
		DistMx[i].resize(SeedCount, FLT_MAX);
		DistMx[i][i] = 0;
		}

	uint PairCount = (SeedCount*(SeedCount-1))/2;
	uint i = 1;
	uint j = 0;
	uint PairCounter = 0;
	uint ThreadCount = GetRequestedThreadCount();
#pragma omp parallel for num_threads(ThreadCount)
	for (int PairIndex = 0; PairIndex < (int) PairCount; ++PairIndex)
		{
		LOCK();
		ProgressStep(PairCounter++, PairCount, "Seed pairs");
		UNLOCK();
		uint SeqIndexi = m_SeedSeqIndexes[i];
		uint SeqIndexj = m_SeedSeqIndexes[j];
		double d = GetProtDist(SeqIndexi, SeqIndexj);
		DistMx[i][j] = float(d);
		DistMx[j][i] = float(d);

		++j;
		assert(j <= i);
		if (j == i)
			{
			++i;
			j = 0;
			}
		}

	for (uint i = 0; i < SeedCount; ++i)
		{
		uint SeqIndex = m_SeedSeqIndexes[i];
		const char *SeqLabel = 
		  m_Builder->m_Seqs->GetLabel(SeqIndex);
		string Label;
		Ps(Label, "Seed%u.%s", i, SeqLabel);
		m_UPGMA.m_Labels.push_back(Label);
		//Log("%s\n", Label.c_str());
		}

	m_UPGMA.Run(LINKAGE_Avg, m_Tree);

	m_Children.clear();
	m_Children.resize(SeedCount);
	for (uint i = 0; i < SeedCount; ++i)
		{
		uint SeedSeqIndex = m_SeedSeqIndexes[i];
		m_Children[i] = new GTBNode;
		m_Children[i]->m_Builder = m_Builder;
		m_Children[i]->m_SeqIndexes.push_back(SeedSeqIndex);
		}

	uint Counter = 0;
#pragma omp parallel for num_threads(ThreadCount)
	for (int i = 0; i < (int) SeqCount; ++i)
		{
		LOCK();
		ProgressStep(Counter++, SeqCount, "Assigning");
		UNLOCK();
		uint SeqIndex = m_SeqIndexes[i];
		if (SeedSet.find(SeqIndex) != SeedSet.end())
			continue;
		double Mind = DBL_MAX;
		uint BestSeedIndex = UINT_MAX;
		for (uint SeedIndex = 0; SeedIndex < SeedCount; ++SeedIndex)
			{
			uint SeedSeqIndex = m_SeedSeqIndexes[SeedIndex];
			double d = GetProtDist(SeqIndex, SeedSeqIndex);
			if (d < Mind)
				{
				Mind = d;
				BestSeedIndex = SeedIndex;
				}
			}
		Log("Assign %s ==> %s\n",
		  m_Builder->m_Seqs->GetLabel(SeqIndex),
		  m_UPGMA.m_Labels[BestSeedIndex].c_str());
		asserta(BestSeedIndex < SIZE(m_Children));
		m_Children[BestSeedIndex]->m_SeqIndexes.push_back(SeqIndex);
		}

	for (uint ChildIndex = 0; ChildIndex < SeedCount; ++ChildIndex)
		m_Children[ChildIndex]->Run();
	}

#endif // 0
