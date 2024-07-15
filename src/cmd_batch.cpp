#include "muscle.h"
#include "bench.h"
#include "m3alnparams.h"

void cmd_batch()
	{
	M3AlnParams AP;
	AP.SetFromCmdLine(false);

	vector<string> Names;
	ReadStringsFromFile(opt(batch), Names);
	string InputDir = opt(indir);
	string OutputDir = opt(outdir);
	Dirize(InputDir);
	Dirize(OutputDir);
	opt_quiet = true;
	optset_quiet = true;

	const uint ThreadCount = GetRequestedThreadCount();
	vector<Muscle3 *> M3s;
	for (uint i = 0; i < ThreadCount; ++i)
		M3s.push_back(new Muscle3);

	const uint N = SIZE(Names);
	uint Counter = 0;
#pragma omp parallel for num_threads(ThreadCount)
	for (int i = 0; i < (int) N; ++i)
		{
		const string &Name = Names[i];
		const string &FastaFileName = InputDir + Name;
		uint ThreadIndex = GetThreadIndex();
		Muscle3 *M3 = M3s[ThreadIndex];

#pragma omp critical
		{
		opt_quiet = false;
		ProgressStep(Counter++, N, "Aligning %u sets %s",
		  N, FastaFileName.c_str());
		opt_quiet = true;
		}

		MultiSequence *MS = new MultiSequence;
		MS->LoadMFA(FastaFileName, true);

		string OutputFileName = OutputDir + Name;

		M3->Run(AP, *MS);
		M3->m_FinalMSA->WriteMFA(OutputFileName);

		delete MS;
		}
	}
