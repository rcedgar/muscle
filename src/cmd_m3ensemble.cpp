#include "muscle.h"
#include "bench.h"
#include "m3alnparams.h"
#include "muscle3.h"

// -replicates 4  testdir=d:/int/balibase/aln/mawc_m3efa/ n=386 avgq=0.8385       avgtc=0.5463
// -replicates 16 testdir=d:/int/balibase/aln/mawc_m3efa/ n=386 avgq=0.8450       avgtc=0.5500
// -replicates 64 testdir=d:/int/balibase/aln/mawc_m3efa/ n=386 avgq=0.8466       avgtc=0.5548

void cmd_m3ensemble()
	{
	asserta(!optset_n);
	const string &FastaFN = g_Arg1;
	asserta(optset_output);
	const string &OutputFN = opt(output);

	FILE *fOut = CreateStdioFile(OutputFN);

	MultiSequence *MS = new MultiSequence;
	MS->LoadMFA(FastaFN, true);

	opt_quiet = true;
	optset_quiet = true;

	vector<Muscle3 *> M3s;
	vector<M3AlnParams *> APs;
	const float Delta = 0.1f;
	uint ThreadCount = GetRequestedThreadCount();
	for (uint i = 0; i < ThreadCount; ++i)
		{
		Muscle3 *M3 = new Muscle3;
		M3AlnParams *AP = new M3AlnParams;

		M3s.push_back(M3);
		APs.push_back(AP);
		}
	int Replicates = 16;
	if (optset_replicates)
		Replicates = int(opt(replicates));
	asserta(Replicates > 0);

#pragma omp parallel for num_threads(ThreadCount)
	for (int ii = 0; ii < Replicates; ++ii)
		{
		uint ThreadIndex = GetThreadIndex();
		Muscle3 &M3 = *M3s[ThreadIndex];
		M3AlnParams &AP = *APs[ThreadIndex];

		uint i = uint(ii);
		uint PerturbSeed = i/4;
		uint ParamGroup = (Replicates == 4 ? 0 : (i*7)%4);
		uint PctId = UINT_MAX;
		switch (i%4)
			{
		case 0: PctId = 90; break;
		case 1: PctId = 80; break;
		case 2: PctId = 70; break;
		case 3: PctId = 62; break;
		default: asserta(false);
			}

		AP.SetBlosum(PctId, ParamGroup, FLT_MAX, FLT_MAX, PerturbSeed,
		  Delta, Delta, Delta, true);
		M3.Run(AP, *MS);
#pragma omp critical
		{
		fprintf(fOut, "<blosum%u:%u.perturb%u.delta%.3g\n",
		  PctId, ParamGroup, PerturbSeed, Delta);
		M3.m_FinalMSA->WriteMFA(fOut);
		}
		}

	CloseStdioFile(fOut);
	}
