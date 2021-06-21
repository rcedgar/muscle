#include "myutils.h"
#include "muscle.h"

void LogTransProbs(const PairHMM &HMM)
	{
	float SumInitSingle = 0;
	for (uint i = 0; i < 20; ++i)
		SumInitSingle += emitSingleDefault[i];
	asserta(feq(SumInitSingle, 1.0f));

	float SumInitPairs = 0;
	for (uint i = 0; i < 20; ++i)
		for (uint j = 0; j < 20; ++j)
			{
			float Pij = emitPairsDefault[i][j];
			float Pji = emitPairsDefault[j][i];
			if (Pij == 0)
				{
				asserta(Pji != 0);
				SumInitPairs += Pji;
				}
			else
				{
				asserta(Pij != 0);
				SumInitPairs += Pij;
				}
			}
	asserta(feq(SumInitPairs, 1.0f));

// I1 transitions
	float Trans_M_to_I1 =  gapOpen2Default[0];
	float Trans_I1_to_I1 =  gapExtend2Default[0];
	float Trans_I1_to_M = 1 - Trans_I1_to_I1;

// I2 transitions
	float Trans_M_to_I2 =  gapOpen2Default[2];
	float Trans_I2_to_I2 =  gapExtend2Default[2];
	float Trans_I2_to_M = 1 - Trans_I2_to_I2;

	float Trans_M_to_M = 1 - 2*Trans_M_to_I1 - 2*Trans_M_to_I2;

	ProgressLog("\n");

	ProgressLog("%7.5f  M -> M\n", Trans_M_to_M);
	ProgressLog("%7.5f  M -> I1\n", Trans_M_to_I1);
	ProgressLog("%7.5f  M -> I2\n", Trans_M_to_I2);

	ProgressLog("\n");
	ProgressLog("%7.5f  I1 -> M\n", Trans_I1_to_M);
	ProgressLog("%7.5f  I1 -> I1\n", Trans_I1_to_I1);

	ProgressLog("\n");
	ProgressLog("%7.5f  I2 -> M\n", Trans_I2_to_M);
	ProgressLog("%7.5f  I2 -> I2\n", Trans_I2_to_I2);

	ProgressLog("\n");
	ProgressLog("transProb:\n");
	for (uint i = 0; i < uint(NumMatrixTypes); ++i)
		{
		ProgressLog("[%u]  ", i);
		float Sum = 0;
		for (uint j = 0; j < uint(NumMatrixTypes); ++j)
			{
			float P = EXP(HMM.transProb[i][j]);
			ProgressLog("  %7.5f", P);
			Sum += P;
			}
		ProgressLog("  (%.5f)\n", Sum);
		}
	ProgressLog("\n");
	}

static void GetTweaks(float &dM_I1, float &dI1_I1, 
  float &dM_I2, float &dI2_I2)
	{
	uint PctA = randu32()%101;
	uint PctB = randu32()%101;
	uint PctC = randu32()%101;
	uint PctD = randu32()%101;

	bool BoostA = (PctA%2 == 0);
	bool BoostB = (PctB%2 == 0);
	bool BoostC = (PctC%2 == 0);
	bool BoostD = (PctD%2 == 0);

	float fA = (BoostA ? float(200)/(PctA + 100) : float(PctA)/100);
	float fB = (BoostB ? float(200)/(PctB + 100) : float(PctB)/100);
	float fC = (BoostC ? float(200)/(PctC + 100) : float(PctC)/100);
	float fD = (BoostD ? float(200)/(PctD + 100) : float(PctD)/100);

	dM_I1	= fA;
	dI1_I1	= fB;
	dM_I2	= fC;
	dI2_I2	= fD;
	}

void cmd_explore()
	{
	const string &InputFileName = opt(explore);
	MultiSequence InputSeqs;
	InputSeqs.FromFASTA(InputFileName);

	const PairHMM &HMM = InitProbcons();
	LogTransProbs(HMM);

	const uint ITERS = 10;
	for (uint Iter = 0; Iter < ITERS; ++Iter)
		{
		ProgressLog("\n");
		ProgressLog("__________________ Iter %u __________\n", Iter);

		PairHMM HMM;
		HMM.SetProbs2(0.815e-5f, 0.159f, 0.0119f, 0.008f, 0.397f, 0.899f);
		LogTransProbs(HMM);
		}
	}
