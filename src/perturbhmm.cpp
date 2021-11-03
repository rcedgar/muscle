#include "muscle.h"

static void Perturb(float &P, float Var)
	{
	asserta(Var >= 0 && Var < 1);
	uint Pct = randu32()%100;
	float Fract = Pct/100.0f;
	asserta(Fract >= 0 && Fract <= 1);
	float Lo = 1.0f - Var;
	float Hi = 1.0f + Var;
	float d = Lo + (Hi - Lo)*Fract;
	P *= d;
	}

void HMMParams::PerturbProbs(uint Seed)
	{
	if (Seed == 0)
		return;

	ResetRand(Seed);
	asserta(m_Var > 0 && m_Var < 1);

	for (uint i = 0; i < SIZE(m_Trans); ++i)
		Perturb(m_Trans[i], m_Var);

	const uint AlphaSize = GetAlphaSize();
	for (uint i = 0; i < AlphaSize; ++i)
		for (uint j = 0; j <= i; ++j)
			{
			float P = m_Emits[i][j];
			Perturb(P, m_Var);
			m_Emits[i][j] = P;
			m_Emits[j][i] = P;
			}
	Normalize();
	}

void HMMParams::Compare(const HMMParams &HP1, const HMMParams &HP2,
  float &MeanTransDelta, float &MeanEmitDelta)
	{
	const uint NT = SIZE(HP1.m_Trans);
	const uint AlphaSize = HP1.GetAlphaSize();
	asserta(SIZE(HP2.m_Trans) == NT);
	asserta(SIZE(HP1.m_Emits) == AlphaSize);
	asserta(SIZE(HP2.m_Emits) == AlphaSize);
	float SumT = 0;
	for (uint i = 0; i < NT; ++i)
		{
		float P1 = HP1.m_Trans[i];
		float P2 = HP2.m_Trans[i];
		SumT += abs(P1 - P2);
		}
	MeanTransDelta = SumT/NT;

	float SumE = 0;
	for (uint i = 0; i < AlphaSize; ++i)
		{
		for (uint j = 0; j < AlphaSize; ++j)
			{
			float P1 = HP1.m_Emits[i][j];
			float P2 = HP2.m_Emits[i][j];
			SumE += abs(P1 - P2);
			}
		}
	MeanEmitDelta = SumE/(AlphaSize*AlphaSize);
	}

static void Run1(uint Iter)
	{
	bool Nucleo = opt(nt);
	HMMParams HPDef;
	HPDef.FromDefaults(Nucleo);

	HMMParams HP;
	HP.FromDefaults(Nucleo);
	HP.PerturbProbs(Iter);

	float MeanTransDelta;
	float MeanEmitDelta;
	HMMParams::Compare(HPDef, HP, MeanTransDelta, MeanEmitDelta);

	ProgressLog("Iter %u, trans %8.6f, emit %8.6f\n",
	  Iter, MeanTransDelta, MeanEmitDelta);
	}

void cmd_perturbhmm()
	{
	const string &sIters = opt(perturbhmm);
	const uint ITERS = StrToUint(sIters);

	for (uint Iter = 0; Iter < ITERS; ++Iter)
		Run1(Iter);
	}
