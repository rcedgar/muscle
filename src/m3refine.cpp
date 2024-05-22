#include "muscle.h"
#include "muscle3.h"

static void NormalizeWeights(vector<float> &SeqWeights)
	{
	const uint N = SIZE(SeqWeights);
	asserta(N > 0);
	float Sum = 0;
	for (uint i = 0; i < N; ++i)
		Sum += SeqWeights[i];
	for (uint i = 0; i < N; ++i)
		SeqWeights[i] /= Sum;
	}

static void SplitIndexes3(uint N,
  vector<vector<uint> > &IndexVec)
	{
	IndexVec.clear();
	IndexVec.resize(3);

	uint Ix0 = randu32()%(N-1);
	uint Ix1 = randu32()%(N-1);
	if (Ix1 == Ix0)
		Ix1 = (Ix1 + 1)%(N-1);
	if (Ix0 > Ix1)
		swap(Ix0, Ix1);

	for (uint i = 0; i <= Ix0; ++i)
		IndexVec[0].push_back(i);

	for (uint i = Ix0 + 1; i <= Ix1; ++i)
		IndexVec[1].push_back(i);

	for (uint i = Ix1 + 1; i < N; ++i)
		IndexVec[2].push_back(i);
#if 0
	{
	for (uint i = 0; i < 3; ++i)
		{
		const uint n = SIZE(LabelVec[i]);
		Log(" [%u]=", i);
		for (uint j = 0; j < n; ++j)
			Log(" %u", LabelVec[i][j]);
		}
	Log("\n");
	}
#endif
	}

void M3Refine(const MultiSequence &InputMSA,
  const M3AlnParams &AP, const vector<float> &SeqWeights,
  MultiSequence &RefinedMSA)
	{
	const Mx2020 &SubstMx_Letter = AP.m_SubstMx_Letter;
	const float GapOpen = AP.m_GapOpen;

	const uint SeqCount = InputMSA.GetSeqCount();
	MultiSequence CurrentBestMSA;
	CurrentBestMSA.Copy(InputMSA);

	Profile3 Prof3;
	Prof3.FromMSA(CurrentBestMSA, SubstMx_Letter, GapOpen, SeqWeights);
	float CurrentBestSelfScore = Prof3.GetSelfScore();

	CacheMem3 CM;
	for (uint Iter = 0; Iter < 32; ++Iter)
		{
		vector<vector<uint> > IndexVec;
		SplitIndexes3(SeqCount, IndexVec);

		MultiSequence *SubMSA0 = CurrentBestMSA.Project(IndexVec[0]);
		MultiSequence *SubMSA1 = CurrentBestMSA.Project(IndexVec[1]);
		MultiSequence *SubMSA2 = CurrentBestMSA.Project(IndexVec[2]);

		vector<float> SeqWeights0;
		vector<float> SeqWeights1;
		vector<float> SeqWeights2;

		for (uint i = 0; i < SIZE(IndexVec[0]); ++i)
			{
			uint SeqIndex = IndexVec[0][i];
			float Weight = SeqWeights[SeqIndex];
			SeqWeights0.push_back(Weight);
			}
		for (uint i = 0; i < SIZE(IndexVec[1]); ++i)
			{
			uint SeqIndex = IndexVec[1][i];
			float Weight = SeqWeights[SeqIndex];
			SeqWeights1.push_back(Weight);
			}
		for (uint i = 0; i < SIZE(IndexVec[2]); ++i)
			{
			uint SeqIndex = IndexVec[2][i];
			float Weight = SeqWeights[SeqIndex];
			SeqWeights2.push_back(Weight);
			}

		NormalizeWeights(SeqWeights0);
		NormalizeWeights(SeqWeights1);
		NormalizeWeights(SeqWeights2);

		Profile3 Prof0;
		Profile3 Prof1;
		Profile3 Prof2;
		Prof0.FromMSA(*SubMSA0, SubstMx_Letter, GapOpen, SeqWeights0);
		Prof1.FromMSA(*SubMSA1, SubstMx_Letter, GapOpen, SeqWeights1);
		Prof2.FromMSA(*SubMSA2, SubstMx_Letter, GapOpen, SeqWeights2);

		string Path01;
		string Path02;
		string Path12;
		NWSmall3(CM, Prof0, Prof1, Path01);
		NWSmall3(CM, Prof0, Prof2, Path02);
		NWSmall3(CM, Prof1, Prof2, Path12);

		Log("\n");
		Log("Path01=%s\n", Path01.c_str());
		Log("Path02=%s\n", Path01.c_str());
		Log("Path12=%s\n", Path01.c_str());

		delete SubMSA0;
		delete SubMSA1;
		delete SubMSA2;
		}
	}

void cmd_m3refine()
	{
	MultiSequence MSA;
	MSA.FromFASTA(g_Arg1);
	asserta(MSA.IsAligned());
	const uint SeqCount = MSA.GetSeqCount();

	vector<string> Labels;
	vector<float> SeqWeights;
	float w = 1.0f/SeqCount;
	for (uint i = 0; i < SeqCount; ++i)
		{
		string Label = (string) MSA.GetLabel(i);
		Labels.push_back(Label);
		SeqWeights.push_back(w);
		}

	vector<vector<float> > DistMx;
	GetKimuraDistMx(MSA, DistMx);

	UPGMA5 U5;
	Tree T;
	U5.Init(Labels, DistMx);
	U5.Run(LINKAGE_Biased, T);

	ClustalWeights CW;
	CW.Run(MSA, T, SeqWeights);

	M3AlnParams AP;
	AP.SetFromCmdLine(false);

	MultiSequence RefinedMSA;
	M3Refine(MSA, AP, SeqWeights, RefinedMSA);
	}
