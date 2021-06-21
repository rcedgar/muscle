#include "muscle.h"

static float GetRandDelta_Trans()
	{
	uint Pct = randu32()%100;
	float Fract = Pct/100.0f;
	float r = 1.6f*Fract;
	float d = 0.2f + r;
	return d;
	}

static float GetRandDelta_Emit(float Var)
	{
	asserta(Var >= 0 && Var < 1);
	uint Pct = randu32()%100;
	float Fract = Pct/100.0f;
	asserta(Fract >= 0 && Fract <= 1);
	float Lo = 1.0f - Var;
	float Hi = 1.0f + Var;
	float d = Lo + (Hi - Lo)*Fract;
	return d;
	}

void HMMParams::Delta()
	{
	float dStartIS = GetRandDelta_Trans();
	float dStartIL = GetRandDelta_Trans();
	float dShortOpen = GetRandDelta_Trans();
	float dShortExtend = GetRandDelta_Trans();
	float dLongOpen = GetRandDelta_Trans();
	float dLongExtend = GetRandDelta_Trans();
	DeltaTransProbs(dStartIS, dStartIL, dShortOpen, dShortExtend,
	  dLongOpen, dLongExtend);

	vector<vector<float> > Deltas(HMM_ALPHASIZE);
	for (uint i = 0; i < HMM_ALPHASIZE; ++i)
		{
		char c = HMM_ALPHA[i];
		for (uint j = 0; j < HMM_ALPHASIZE; ++j)
			{
			float d = GetRandDelta_Emit(0.3f);
			Deltas[i].push_back(d);
			}
		}
	DeltaEmitProbs(Deltas);
	}

static void Run1(MultiSequence &InputSeqs, uint Iter, const string &OutDir)
	{
	HMMParams HP;
	HP.FromDefaults();
	HP.Delta();
	HP.ToPairHMM();

	string ParamsFileName;
	Ps(ParamsFileName, "%s/hmm%u.tsv", OutDir.c_str(), Iter);
	HP.ToFile(ParamsFileName);

	string MSAFileName;
	Ps(MSAFileName, "%s/msa%u.afa", OutDir.c_str(), Iter);
	MultiSequence *MSA = RunMPC(&InputSeqs);

	MSA->WriteMFA(MSAFileName);
	delete MSA;
	}

void cmd_deltahmm()
	{
	const string &InputFileName = opt(deltahmm);
	string OutDir = opt(outdir);
	Dirize(OutDir);

	MultiSequence InputSeqs;
	InputSeqs.FromFASTA(InputFileName, true);

	InitProbcons();

	const uint ITERS = 10;
	for (uint Iter = 0; Iter < ITERS; ++Iter)
		{
		Progress("\n");
		Progress("__________________ Iter %u / %u _______________\n",
		  Iter + 1, ITERS);
		Progress("\n");
		Run1(InputSeqs, Iter, OutDir);
		}
	}
