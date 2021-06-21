#include "muscle.h"
#include "hmmparams.h"

static bool g_InitDone = false;

void InitProbcons()
	{
	if (g_InitDone)
		return;

	SetAlpha(ALPHA_Amino);

	HMMParams HP;
	if (optset_hmmin)
		{
		const string FileName = opt(hmmin);
		ProgressLog("Reading HMM parameters from %s\n", FileName.c_str());
		HP.FromFile(FileName);
		}
	else
		HP.FromDefaults();

	if (optset_perturb)
		{
		uint Seed = opt(perturb);
		if (Seed > 0)
			{
			ProgressLog("Perturbing HMM parameters with seed %us\n", Seed);
			ResetRand(Seed);
			HP.Delta();
			}
		}

	HP.ToPairHMM();
	if (optset_hmmout)
		HP.ToFile(opt(hmmout));

	HP.ToPairHMM();

	g_InitDone = true;
	}
