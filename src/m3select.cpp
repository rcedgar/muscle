#include "muscle.h"
#include "bench.h"
#include "m3alnparams.h"
#include "muscle3.h"

/***
Default Muscle3 progressive
testdir=d:/int/balibase/aln/m3/ n=386 avgq=0.8354       avgtc=0.5412

m3select variants
[bafeef2] testdir=d:/int/balibase/aln/m3select/ n=386 avgq=0.8418 avgtc=0.5444 16 replicates, various substmxs, perturb everything
[613b055] testdir=d:/int/balibase/aln/m3select/ n=386 avgq=0.8397 avgtc=0.5456 16 replicates, B62:0 only, perturb distmx only
[859cadb] testdir=d:/int/balibase/aln/m3select/ n=386 avgq=0.8354 avgtc=0.5381 64 replicates, B62:0 only, perturb distmx only
***/

void cmd_m3select()
	{
	const string &FastaFN = g_Arg1;

	MultiSequence *MS = new MultiSequence;
	MS->LoadMFA(FastaFN, true);

	opt_quiet = true;
	optset_quiet = true;

	M3AlnParams MasterAP;
	MasterAP.SetBlosum(62, 0, FLT_MAX, FLT_MAX, 0, 0, 0, 0, false);
	const Mx2020 &Master_SubstMx_Letter = MasterAP.m_SubstMx_Letter;
	const float Master_GapOpen = MasterAP.m_GapOpen;

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
	int Replicates = 64;
	if (optset_replicates)
		Replicates = int(opt(replicates));
	asserta(Replicates > 0);

	const MultiSequence *BestMSA = 0;
	float BestSelfScore = 0;

#pragma omp parallel for num_threads(ThreadCount)
	for (int ii = 0; ii < Replicates; ++ii)
		{
		uint ThreadIndex = GetThreadIndex();
		Muscle3 &M3 = *M3s[ThreadIndex];
		M3AlnParams &AP = *APs[ThreadIndex];

		uint PerturbSeed = uint(ii);
		uint ParamGroup = 0;
		uint PctId = 62;

		float PerturbSubstMxDelta = 0;
		float PerturbGapParamsDelta = 0;
		float PerturbDistMxDelta = 0.1f;

		AP.SetBlosum(PctId, ParamGroup, FLT_MAX, FLT_MAX, PerturbSeed,
		  PerturbSubstMxDelta, PerturbGapParamsDelta, PerturbDistMxDelta, false);
		M3.Run(AP, *MS);
		Profile3 Prof;
		Prof.FromMSA(*M3.m_FinalMSA,
		  Master_SubstMx_Letter, Master_GapOpen,
		  M3.m_InputSeqWeights);
		float SelfScore = Prof.GetSelfScore();
#pragma omp critical
		{
		if (SelfScore > BestSelfScore)
			{
			BestMSA = M3.m_FinalMSA;
			BestSelfScore = SelfScore;
			}
		}
		}
	asserta(BestMSA != 0);
	BestMSA->ToFasta(opt(output));
	}
