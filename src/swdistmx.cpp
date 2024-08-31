#include "muscle.h"
#include "mx.h"
#include "tracebit.h"
#include "xdpmem.h"
#include "swtrace.h"
#include "upgma5.h"
#include <mutex>

float SWFast_Seqs_BLOSUM62(XDPMem &Mem, const Sequence &A, const Sequence &B,
  float Open, float Ext, uint &Loi, uint &Loj, uint &Leni, uint &Lenj,
  string &Path);

struct ThreadData
	{
	const MultiSequence *Input = 0;
	Mx<float> *DistMx = 0;
	uint PairIndex = UINT_MAX;
	uint Idx1 = UINT_MAX;
	uint Idx2 = UINT_MAX;
	float Open = FLT_MAX;
	float Ext = FLT_MAX;
	mutex Lock;
	};

static void ThreadBody(uint ThreadIndex, void *ptrData)
	{
	ThreadData &TD = *(ThreadData *) ptrData;
	Mx<float> &DistMx = *TD.DistMx;
	const MultiSequence &Input = *TD.Input;
	const uint SeqCount = Input.GetSeqCount();
	uint PairCount = (SeqCount*(SeqCount - 1))/2;
	const float Open = TD.Open;
	const float Ext = TD.Ext;
	XDPMem Mem;
	for (;;)
		{
		TD.Lock.lock();
		bool Done = (TD.PairIndex == PairCount);
		uint Idx1 = TD.Idx1;
		uint Idx2 = TD.Idx2;
		if (!Done)
			{
			ProgressStep(TD.PairIndex++, PairCount, "Aligning");
			asserta(Idx1 < SeqCount);
			asserta(Idx2 < SeqCount);
			asserta(Idx1 != Idx2);
			++TD.Idx2;
			if (TD.Idx2 == SeqCount)
				{
				++TD.Idx1;
				if (TD.Idx1 != SeqCount)
					TD.Idx2 = TD.Idx1 + 1;
				}
			}
		TD.Lock.unlock();
		if (Done)
			return;
		asserta(Idx1 < SeqCount);
		asserta(Idx2 < SeqCount);

		const Sequence &Seqi = *Input.GetSequence(Idx1);
		const Sequence &Seqj = *Input.GetSequence(Idx2);

		const uint Li = Seqi.GetLength();
		const uint Lj = Seqj.GetLength();

		const string &Labeli = Seqi.m_Label;
		const string &Labelj = Seqj.m_Label;

		string Path;
		uint Loi, Loj, Leni, Lenj;
		uint ThreadIndex = GetThreadIndex();
		float SWScore = SWFast_Seqs_BLOSUM62(Mem, Seqi, Seqj, Open, Ext,
			Loi, Loj, Leni, Lenj, Path);
		float MeanLength = (Li + Lj)/2.0f;
		float NormScore = SWScore/MeanLength;

		TD.Lock.lock();
		DistMx.Put(Idx1, Idx2, NormScore);
		DistMx.Put(Idx2, Idx1, NormScore);
		Log("%10.10s  %10.10s  %10.1f  %7.4f  %s\n",
			Labeli.c_str(), Labelj.c_str(), SWScore, NormScore, Path.c_str());
		TD.Lock.unlock();
		}
	}

void CalcGuideTree_SW_BLOSUM62(const MultiSequence &Input, Tree &T)
	{
	const uint SeqCount = Input.GetSeqCount();
	vector<string> Labels;
	for (uint i = 0; i < SeqCount; ++i)
		{
		const string &Label = Input.GetLabelStr(i);
		Labels.push_back(Label);
		}

	Mx<float> DistMx;
	DistMx.Alloc("SWDistMx", SeqCount, SeqCount);
	DistMx.PutAll(FLT_MAX);

	ThreadData TD;
	TD.DistMx = &DistMx;
	TD.Input = &Input;
	TD.Idx1 = 0;
	TD.Idx2 = 1;
	TD.Open = -11;
	TD.Ext = -1;
	TD.PairIndex = 0;

	RunThreads(ThreadBody, &TD);

	asserta(TD.Idx1 == SeqCount - 1 && TD.Idx2 == SeqCount);
	Progress("Done.\n");

	vector<vector<float> > DistMxv(SeqCount);
	for (uint i = 0; i < SeqCount; ++i)
		DistMxv[i].resize(SeqCount);

	for (uint i = 0; i < SeqCount; ++i)
		for (uint j = 0; j < SeqCount; ++j)
			DistMxv[i][j] = DistMx.Get(i, j);

	UPGMA5 U;
	U.Init(Labels, DistMxv);
	U.ScaleDistMx();
	U.Run(LINKAGE_Avg, T);
	}

void cmd_swdistmx()
	{
	MultiSequence Input;
	Input.FromFASTA(g_Arg1, true);
	Tree T;
	CalcGuideTree_SW_BLOSUM62(Input, T);
	T.ToFile(opt(guidetreeout));
	}
