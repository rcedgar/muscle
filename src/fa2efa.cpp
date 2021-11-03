#include "muscle.h"
#include "ensemble.h"

void cmd_fa2efa()
	{
	const string &InputFileName = opt(fa2efa);
	const string &OutputFileName = opt(output);

	Ensemble E;
	E.FromMSAPaths(InputFileName);

	Progress("%u seqs, %u MSAs\n", E.GetSeqCount(), E.GetMSACount());
	Progress("Writing %s ...\n", OutputFileName.c_str());
	E.ToEFA(OutputFileName);
	Progress("done.\n");
	}
