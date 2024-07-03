#include "muscle.h"
#include "hmmparams.h"

static bool g_InitDone = false;

void InitProbcons()
	{
	if (g_InitDone)
		return;

	asserta(g_Alpha == ALPHA_Amino || g_Alpha == ALPHA_Nucleo);

	HMMParams HP;
	if (optset_hmmin)
		{
		const string FileName = opt(hmmin);
		ProgressLog("Reading HMM parameters from %s\n", FileName.c_str());
		HP.FromFile(FileName);
		}
	else
		{
		bool Nucleo = (g_Alpha == ALPHA_Nucleo);
		HP.FromDefaults(Nucleo);
		}

	if (optset_perturb)
		{
		uint Seed = opt(perturb);
		if (Seed > 0)
			{
			ProgressLog("Perturbing HMM parameters with seed %u\n", Seed);
			ResetRand(Seed);
			HP.PerturbProbs(Seed);
			}
		}

	if (optset_hmmout)
		HP.ToFile(opt(hmmout));

	HP.CmdLineUpdate();
	HP.ToPairHMM();

	g_InitDone = true;
	}
