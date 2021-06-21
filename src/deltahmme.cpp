#include "muscle.h"

static float GetRandDelta(float Var)
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

static void Run1(MultiSequence &InputSeqs, uint Iter, const string &OutDir)
	{
	HMMParams HP;
	HP.FromDefaults();

	vector<vector<float> > Deltas(HMM_ALPHASIZE);
	for (uint i = 0; i < HMM_ALPHASIZE; ++i)
		{
		char c = HMM_ALPHA[i];
		for (uint j = 0; j < HMM_ALPHASIZE; ++j)
			{
			float d = GetRandDelta(0.3f);
			Deltas[i].push_back(d);
			}
		}
	HP.DeltaEmitProbs(Deltas);

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

void cmd_deltahmme()
	{
	const string &InputFileName = opt(deltahmme);
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
