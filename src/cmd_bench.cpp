#include "muscle.h"
#include "bench.h"
#include "m3alnparams.h"

void cmd_bench()
	{
	M3AlnParams AP;
	AP.SetFromCmdLine(false);

	Bench B;
	string RefDir = opt(refdir);
	Dirize(RefDir);
	B.Load(g_Arg1, RefDir);
	B.m_ShowProgress = true;
	uint MSACount = SIZE(B.m_Inputs);

	optset_quiet = true;
	opt_quiet = true;
	B.Run(AP);
	B.TCsToFile(opt(tsvout));

	optset_quiet = false;
	opt_quiet = false;
	ProgressLog("AvgQ=%.3f AvgTC=%.3f N=%u\n",
		B.m_MeanQ, B.m_MeanTC, MSACount);
	}

void cmd_bench_blosums()
	{
	Bench B;
	string RefDir = opt(refdir);
	Dirize(RefDir);
	B.Load(g_Arg1, RefDir);
	B.m_ShowProgress = true;
	uint MSACount = SIZE(B.m_Inputs);
	FILE *fTsv = CreateStdioFile(opt(tsvout));

	optset_quiet = true;
	opt_quiet = true;

	if (fTsv != 0)
		{
		fprintf(fTsv, "BLOSUM");
		fprintf(fTsv, "\tParamSet");
		fprintf(fTsv, "\tQ");
		fprintf(fTsv, "\tTC");
		fprintf(fTsv, "\tPerturbSeed");
		fprintf(fTsv, "\tDelta");
		fprintf(fTsv, "\n");
		}
	for (uint PerturbSeed = 0; PerturbSeed < 6; ++PerturbSeed)
		{
		float Delta = 0.05f*PerturbSeed;

		for (uint k = 0; k < 4; ++k)
			{
			uint PctId = UINT_MAX;
			switch (k)
				{
			case 0: PctId = 90; break;
			case 1: PctId = 80; break;
			case 2: PctId = 70; break;
			case 3: PctId = 62; break;
				}

			for (uint n = 0; n < 4; ++n)
				{
				M3AlnParams AP;
				AP.SetBlosum(PctId, n, FLT_MAX, FLT_MAX,
				  PerturbSeed, Delta, Delta, Delta, true);
				B.Run(AP);

				printf("BLOSUM%u:%u perturb=%u delta=%7.3g AvgQ=%.4f AvgTC=%.4f N=%u\n",
					PctId, n, PerturbSeed, Delta, B.m_MeanQ, B.m_MeanTC, MSACount);
				Log("BLOSUM%u:%u  perturb=%u delta=%7.3g AvgQ=%.4f AvgTC=%.4f N=%u\n",
					PctId, n, PerturbSeed, Delta, B.m_MeanQ, B.m_MeanTC, MSACount);

				if (fTsv != 0)
					{
					fprintf(fTsv, "%u", PctId);
					fprintf(fTsv, "\t%u", n);
					fprintf(fTsv, "\t%.5f", B.m_MeanQ);
					fprintf(fTsv, "\t%.5f", B.m_MeanTC);
					fprintf(fTsv, "\t%u", PerturbSeed);
					fprintf(fTsv, "\t%.3g", Delta);
					fprintf(fTsv, "\n");
					}
				}
			}
		}
	CloseStdioFile(fTsv);
	}
